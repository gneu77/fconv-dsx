#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<map>
#include<tr1/unordered_map>
#include<algorithm>
#include<cmath>
#include"structure_GN.h"
#include"files_GN.h"
#include"linalg_GN.hpp"
#include"grid_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"
#include"string_fu_GN.hpp"
#include"sphere_grid_coords.hpp"

#include"atom_GN.h"
#include"molecule_GN.h"
#include"protein_GN.h"
#include"structure_additional_GN.h"
#include"atom_properties_GN.h"
#include"octree_GN.hpp"

#include"nr3.h"
#include"mins.h"
#include"mins_ndim.h"

#include<time.h>

using namespace std;
using namespace TXT;


const float crad = 0.03;   //Radius fuer die Cylinder in der Visualisierung
const float min_sphere_rad = 0.02; //Mindestgroesse fuer die Kugeln
const float max_sphere_rad = 6.0; //Maximalgroesse fuer die Kugeln
const float flex_cav_rad = 7.0; //Radius in dem die Cavities fuer flexibles docking ausgeschnitten werden
const float specific_inter_dst = 5.0;

const string version_string = "0.90";
const string version_date = "31.01.2012";
const string copyright_string =
" Copyright (C) 2009, 2010, 2011, 2012 Gerd Neudert and Gerhard Klebe\n\
  Usage of DSX is free without any limitations.\n\
  There is NO warranty, not even for MERCHANTABILITY or FITNESS\n\
  FOR A PARTICULAR PURPOSE.\n\n";

//! Gewichte fuer die einzelnen Scoring-Terme:
float const weight_limit = 0.00001;
float pair_weight = 1.0;
float torsion_weight = 0.0;
float intra_weight = 0.0;
float sas_weight = 0.0;
float wat_weight = 1.0;


const float sas_lig_lig_spacing = 4.;
const float sas_lig_pro_spacing = 2.;
const float MAX_RAD2 = 3.2;
const float LMAX_RAD2 = 3.5;


//! Hier switchen welche SAS-Variante genutzt werden soll:
#define MODE_RATIO
//#define MODE_DELTA // Ist DEUTLICH schlechter


class TRIPLE;
class SAS_POINT;
class SCORER;

//=================================================================================================
// Deklarationen
//=================================================================================================

class TRIPLE {
    public:
        vec3d<float> a;
        vec3d<float> b;
        float val;
        TRIPLE(vec3d<float> &va,vec3d<float> &vb,float vw) : a(va),b(vb),val(vw) {}
        ~TRIPLE() {}
};

class SAS_POINT {
    public:
        vec3d<float> coord;
        bool is_opposite;
        ATOM* closest_pro;
        float pro_dst;
        ATOM* closest_lig;
        float lig_dst;
        SAS_POINT(vec3d<float> const& p):coord(p),is_opposite(false),closest_pro(0),pro_dst(999999.),closest_lig(0),lig_dst(999999.) {}
        ~SAS_POINT() {}
};


inline int64_t float2int64(float const& d) {
    return static_cast<int64_t>(d<0.?d-0.5:d+0.5);
}


class SCORER {
    private:
        //!--------------------------------------------------------------------------------
        //!Die folgenden Parameter werden sich in absehbarer Zeit nicht aendern. Trotzdem
        //!werden sie erstmal variabel gelassen. Die Methode calc_additional_params()
        //!ist somit obligatorisch und sollte sich auch nicht aendern.
        static int max_real_bin; //groesster bin anhand max_dst
        static int max_tors_bin;
        
        static float max_dst; //maximale Distanz fuer das Scoren
        static float max_lim;
        static float max_square_dst;
        static float target_bin_size;
        static float real_bin_size;
        static float tors_bin_size;
        //!--------------------------------------------------------------------------------
        
        //!--------------------------------------------------------------------------------
        //!Uebergebene Parameter:
        string pro_file; //Name des Protein Files (pdb oder mol2)
        string lig_file; //Name des Ligand Files (mol2)
        string cof_file; //Name des optionalen Kofaktor Files (mol2)
        string wat_file; //Name des optionalen Wasser Files (mol2)
        string met_file; //Name des optionalen Metall Files (mol2)
        string ref_file; //Name des optionalen Referenzstruktur Files (mol2)
        string pot_dir; //Name des Verzeichnis, in dem sie die zu benutzenden Potentiale befinden
        string pro_type; //Dateityp des Protein Files
        string lig_type; //Dateityp des Liganden Files
        
        bool coval_check; //Auf kovalent gebundene Liganden pruefen?
        bool score_pairs;
        bool score_torsions; //Torsionswinkel bewerten?
        bool score_intra; //intramolekulare Clashes scoren?
        bool score_sas; //SAS Term berechnen?
        bool gold_water; //Wassermolekuele aus einem GOLD-Docking beruecksichtigen?
        bool jiggle; //Liganden erst minimieren?
        bool jiggle_first_only; //Nur die erste Pose Minimieren?
        bool calc_rmsd; //RMSD-Werte berechnen?
        bool bron_kerbosch;
        static bool silent_mode; //Nur Errors und Warnings ausgeben?
        bool visualization; //Python Files zur Visualisierung in Pymol rausschreiben?
        bool flex_res; //Flexibles Protein?
        bool as_profile;
        
        int interaction_mode; //Festlegen welche Atom-Atom Paarungen moeglich sind
        int sort_mode; //In welcher Reihenfolge sollen die Ergebnisse ausgegeben werden
        static int verbosity;
        
        float vis_scaling; //Skalierungsfaktor fuer die roten und blauen Kugeln
        float vis_lower; //Atom-Atom Interaktionen mit einem Potential unter dieser Schwelle visualisieren
        float vis_upper; //Atom-Atom Interaktionen mit einem Potential oberhalb dieser Schwelle visualisieren
        float vis_tors; //Torsionen mit einem Wert oberhalb dieser Schwelle werden als bad_torsions visualisiert
        //!--------------------------------------------------------------------------------
        
        static vector<stl_ptr<ATOM> > grid_atoms; //Hier kommen die Proteinatome rein
        
        map<int,map<int,vec3d<float> > > flex_map; //! flex_map[ligandnummer][intern_id] = coords
        
        static tr1::unordered_map<string,int> key_map;
        static tr1::unordered_map<string,int> tors_key_map;
        static tr1::unordered_map<string,int> sas_key_map;
        static tr1::unordered_map<string,int> pro_sas_key_map;
        
        static GRID<2,float> sgrid; //! sgrid[key_map[key]][bin] = score  //Die Paar-Potentiale
        static GRID<2,float> tors_sgrid; //! tors_sgrid[key_map[key]][bin] = score  //Die Torsions-Potentiale
        static GRID<2,float> sas_grid;
        static GRID<2,float> pro_sas_grid;
        
        map<int,float> score_map; //! Ligandnummer(in lig_structure):Score
        map<int,float> tors_score_map; //! Ligandnummer(in lig_structure):Score
        map<int,float> sas_score_map;
        map<int,float> intra_score_map;
        map<int,float> per_atom_score;
        map<int,float> per_contact_score;
        map<int,float> per_tors_score;
        map<int,float> rmsd_map; //! Ligandnummer:RMSD
        
        //!==Teil fuer AS-Profile=====================================
        map<int,map<int,float> > res_score; //! Score pro Residue ( lignum:resnum:score )
        map<int,string> res_name; //! resnum:resname
        //!===========================================================
        
        static float spacing;
        static float hsp;

        int max_ligatm; //Anzahl der Atome im groessten Liganden (nur fuer Visualisierung)
        int max_cofatm;
        int max_watatm;
        int max_metatm;
        
        GRID<2,float> pro_vis; //Scores fuer die Proteinatome
        GRID<2,float> lig_vis; //Scores fuer die Ligandatome
        GRID<2,float> cof_vis; //Scores fuer die Kofaktoratome
        GRID<2,float> wat_vis; //Scores fuer die Wasseratome
        GRID<2,float> met_vis; //Scores fuer die Metallatome
        
        map<int,vector<TRIPLE> > good_con; //Interaktionen mit Werten < vis_lower
        map<int,vector<TRIPLE> > bad_con; //Interaktionen mit Werten > vis_upper
        map<int,vector<int> > bad_tors; //Torsionen mit Werten > vis_tors
        
        static ofstream debug_out;
        static ofstream debug_tors_out;
        static ofstream matrix_out;
        
        //!--------------------------------------------------------------------------------
        //!Structure- und Parser-Objekte fuer die verschiedenen Komponenten:
        stl_ptr<STRUCTURE> pro_structure; //Protein
        stl_ptr<PARSER> pro_parser;
        static stl_ptr<STRUCTURE> lig_structure; //Ligand
        stl_ptr<PARSER> lig_parser;
        stl_ptr<STRUCTURE> cof_structure; //Cofactor
        stl_ptr<PARSER> cof_parser;
        stl_ptr<STRUCTURE> wat_structure; //Wasser
        stl_ptr<PARSER> wat_parser;
        stl_ptr<STRUCTURE> met_structure; //Metalle
        stl_ptr<PARSER> met_parser;
        stl_ptr<STRUCTURE> ref_structure; //Referenz
        stl_ptr<PARSER> ref_parser;
        //!--------------------------------------------------------------------------------
        
        //!--------------------------------------------------------------------------------
        //!Zusaetzliche Attribute fuer die Minimierung:
        static int jig_ligno;
        static vector<vec3d<float> > ref_coords;
        static map<stl_ptr<BOND>,vector<stl_ptr<ATOM> > > r_rot;
        static map<stl_ptr<BOND>,multimap<string,vector<stl_ptr<ATOM> > > > r_tors;
        static vec3d<float> centroid;
        static OCTREE<ATOM>* tree_pro;
        static OCTREE<ATOM>* sas_tree_pro;
        static map<stl_ptr<ATOM>,vector<stl_ptr<ATOM> > > clash_at;
        static int freedom;
        //!--------------------------------------------------------------------------------

        static tr1::unordered_map<string,int> a2h;
        static tr1::unordered_set<string> is_not_reduceable;

        //!--------------------------------------------------------------------------------
        //!Functor zur Berechnung des eigentlichen Scores:
        class ENERGY {
            public:
                double operator()(VecDoub_I &point);
                ENERGY() {}
                ~ENERGY() {}
        };
        //!--------------------------------------------------------------------------------
        
        void clean_mem(); //Reservierten Speicher freigeben
        inline void calc_additional_params();
        inline void load_potentials();
        inline void load_torsion_potentials();
        inline void load_sas_potentials();
        inline void load_pro_sas_potentials();
        void read_gold_water();
        void load_structures();
        void get_r_bonded(stl_ptr<BOND> &to,stl_ptr<ATOM> &curr,tr1::unordered_set<int> &visited);
        inline static bool try_to_replace_with_reduced(string const& A,string const& B,
                                                       string const& C,string const& D,
                                                       string const& a,string const& ckey);
    public:
        SCORER(string &profile,string &ligfile,string &coffile,string &watfile,string &metfile,string &reffile,
               string &potdir,string &protype,string &ligtype,bool covalcheck,bool pairs,bool torsions,bool intra,bool sas,
               bool goldwater,bool jig,bool jigfo,
               bool calcrmsd,bool mcis,bool silentmode,bool vis,bool flexres,bool asp,int interactionmode,
               int sortmode,float visscaling,float vislower,float visupper,float vistors);
        ~SCORER();
        void score();
        void write_results(int n_visi,bool write_txt = true);
};

//!============================================================================================================================
//! Statische Attribute initialisieren:
int SCORER::max_real_bin;
int SCORER::max_tors_bin;
float SCORER::max_dst;
float SCORER::max_lim;
float SCORER::max_square_dst;
float SCORER::target_bin_size;
float SCORER::real_bin_size;
float SCORER::tors_bin_size;
bool SCORER::silent_mode;
int SCORER::verbosity;
vector<stl_ptr<ATOM> > SCORER::grid_atoms;
tr1::unordered_map<string,int> SCORER::key_map;
tr1::unordered_map<string,int> SCORER::tors_key_map;
tr1::unordered_map<string,int> SCORER::sas_key_map;
tr1::unordered_map<string,int> SCORER::pro_sas_key_map;
GRID<2,float> SCORER::sgrid;
GRID<2,float> SCORER::tors_sgrid;
GRID<2,float> SCORER::sas_grid;
GRID<2,float> SCORER::pro_sas_grid;
float SCORER::spacing;
float SCORER::hsp;
ofstream SCORER::debug_out;
ofstream SCORER::debug_tors_out;
ofstream SCORER::matrix_out;
stl_ptr<STRUCTURE> SCORER::lig_structure;
int SCORER::jig_ligno = 0;
vector<vec3d<float> > SCORER::ref_coords;
map<stl_ptr<BOND>,vector<stl_ptr<ATOM> > > SCORER::r_rot;
map<stl_ptr<BOND>,multimap<string,vector<stl_ptr<ATOM> > > > SCORER::r_tors;
vec3d<float> SCORER::centroid;
OCTREE<ATOM>* SCORER::sas_tree_pro;
OCTREE<ATOM>* SCORER::tree_pro;
map<stl_ptr<ATOM>,vector<stl_ptr<ATOM> > > SCORER::clash_at;
int SCORER::freedom;
tr1::unordered_map<string,int> SCORER::a2h;
tr1::unordered_set<string> SCORER::is_not_reduceable;
//!============================================================================================================================

//!============================================================================================================================
SCORER::SCORER(string &profile,string &ligfile,string &coffile,string &watfile,string &metfile,string &reffile,
           string &potdir,string &protype,string &ligtype,bool covalcheck,bool pairs,bool torsions,bool intra,bool sas,
           bool goldwater,bool jig,bool jigfo,
           bool calcrmsd,bool mcis,bool silentmode,bool vis,bool flexres,bool asp,int interactionmode,int sortmode,
           float visscaling,float vislower,float visupper,float vistors) :
           pro_file(profile),lig_file(ligfile),cof_file(coffile),wat_file(watfile),met_file(metfile),ref_file(reffile),
           pot_dir(potdir),pro_type(protype),lig_type(ligtype),coval_check(covalcheck),score_pairs(pairs),score_torsions(torsions),score_intra(intra),
           score_sas(sas),
           gold_water(goldwater),jiggle(jig),jiggle_first_only(jigfo),calc_rmsd(calcrmsd),bron_kerbosch(mcis),visualization(vis),
           flex_res(flexres),as_profile(asp),interaction_mode(interactionmode),sort_mode(sortmode),vis_scaling(visscaling),
           vis_lower(vislower),vis_upper(visupper),vis_tors(vistors) {
    silent_mode = silentmode;
    calc_additional_params();
    sas_tree_pro = 0;
    tree_pro = 0;
    if (score_sas) sphere_grid_coords::initialize();
}

SCORER::~SCORER() {
    clean_mem();
}
//!============================================================================================================================

//!============================================================================================================================
//!Allen reservierten Speicher wieder freigeben:
void SCORER::clean_mem() {    
    if (!pro_structure.zero()) pro_structure.kill();
    if (!pro_parser.zero()) pro_parser.kill();
    if (!lig_structure.zero()) lig_structure.kill();
    if (!lig_parser.zero()) lig_parser.kill();
    if (!wat_structure.zero()) wat_structure.kill();
    if (!wat_parser.zero()) wat_parser.kill();
    if (!met_structure.zero()) met_structure.kill();
    if (!met_parser.zero()) met_parser.kill();
    if (!cof_structure.zero()) cof_structure.kill();
    if (!cof_parser.zero()) cof_parser.kill();
    if (!ref_structure.zero()) ref_structure.kill();
    if (!ref_parser.zero()) ref_parser.kill();
    if (sas_tree_pro) delete sas_tree_pro;
    if (tree_pro) delete tree_pro;
}
//!============================================================================================================================


bool is_metal(string const& s) {
	vector<string> tv;
	string_fu::mysplit(s,tv,'.');
	string lk = tv[0];

	if (lk == "C" || lk == "N" || lk == "O" || lk == "S" || lk == "P" || lk == "H" ||
	    lk == "F" || lk == "Cl" || lk == "Br" || lk == "I" || lk == "B" || lk == "Si" || lk == "X") return false;
	else return true;
}
bool is_hal(string const& s) {
	vector<string> tv;
	string_fu::mysplit(s,tv,'.');
	string lk = tv[0];

	if (lk == "F" || lk == "Cl" || lk == "Br" || lk == "I") return true;
	else return false;
}
string reduce_it(string const& s,tr1::unordered_map<string,int> &hyb_map) {
	if (is_metal(s)) {
        return string("Met");

    } else if (is_hal(s)) {
        return string("Hal");

    } else if (s[0] == 'C') {
		if (hyb_map[s] == 1) return string("C.1");
        else if (hyb_map[s] == 2) return string("C.2");
        else if (hyb_map[s] == 3) return string("C.3");
        else return string("C");

    } else if (s[0] == 'N') {
		if (hyb_map[s] == 1) return string("N.1");
        else if (hyb_map[s] == 2) return string("N.2");
        else if (hyb_map[s] == 3) return string("N.3");
        else return string("N");

    } else if (s[0] == 'O') {
		if (hyb_map[s] == 2) return string("O.2");
        else if (hyb_map[s] == 3) return string("O.3");
        else return string("O");

    } else if (s[0] == 'S') {
        if (s == "Si") return s;
        else return string("S");

    } else if (s[0] == 'P') {
        return string("P");

    } else return string("X");
}

bool SCORER::try_to_replace_with_reduced(string const& A,string const& B,string const& C,
                                         string const& D,string const& a,string const& ckey) {
    if (is_not_reduceable.find(ckey) != is_not_reduceable.end()) return false;

    if (a2h.size() < 1) {
        for (int i=0; i<n_intern_types; ++i) a2h[i_t[i]] = atom_hybridizations[i];
    }

    string Ared = reduce_it(A,a2h);
    if (Ared == "X") return false;
    string Bred = reduce_it(B,a2h);
    if (Bred == "X") return false;
    string Cred = reduce_it(C,a2h);
    if (Cred == "X") return false;
    string Dred = reduce_it(D,a2h);
    if (Dred == "X") return false;

    string red_key = Ared + "_" + Bred + "_" + Cred + "_" + Dred + "_" + a;
    if (Dred < Ared) red_key = Dred + "_" + Cred + "_" + Bred + "_" + Ared + "_" + a;
    else if (Dred == Ared && Cred < Bred) {
        red_key = Dred + "_" + Cred + "_" + Bred + "_" + Ared + "_" + a;
    }

    if (tors_key_map.find(red_key) == tors_key_map.end()) {
        is_not_reduceable.insert(ckey);
        return false;
    }

    tors_key_map[ckey] = tors_key_map[red_key];
    return true;
}

//!============================================================================================================================
//!Obligatorisch (siehe oben)
void SCORER::calc_additional_params() {
    real_bin_size = 0.01;
    max_dst = 6.0;
    max_square_dst = max_dst * max_dst;
    target_bin_size = 0.01;
    
    const float PI = acos(-1.);
    tors_bin_size = PI / 90.; //! 5� Schritte
    max_tors_bin = 89;
    
    max_real_bin = int(max_dst / real_bin_size);
}
//!============================================================================================================================

//!============================================================================================================================
//!Potentiale einladen:
void SCORER::load_potentials() {
    if (!silent_mode) cout << " --> loading pair potentials" << endl;
    string ptmp = pot_dir + "/potentials.keys";
    ifstream f_in;
    string row,key_name,last_key;
    f_in.open(ptmp.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << ptmp << "'" << endl;
        clean_mem();
        exit(1);
    }
    getline(f_in,row);
    f_in.close();
    istringstream is;
    is.str(row);
    int key_num;
    last_key = "*";
    int n_diff = 0;
    while (!is.eof()) {
        key_name = last_key;
        is >> key_name;
        if (key_name == last_key) break; //!wichtig, um nicht fuer den letzten key 2mal einen eintrag zu machen!!!
        key_num = -1;
        is >> key_num;
        
        if (key_num == -1) {
            cerr << c_message<cERROR>(ptmp) << " is not valid (only potentials newer than 08/2009 are valid)" << endl;
            clean_mem();
            exit(1);
        }
        
        key_map[key_name] = key_num;
        if (key_num > n_diff) n_diff = key_num;
    }
    string shelp = "/potentials_repulsive.bin";
    ptmp = pot_dir + shelp;
    sgrid.load(ptmp.c_str(),true);
    max_dst = real_bin_size * (sgrid.get_sizes()[1] - 1); //! max_dst wird auch gebraucht!
    max_square_dst = max_dst * max_dst;
//    if (!silent_mode) {
//        cout << "    -> " << key_map.size() << " different pair potentials loaded from '" << pot_dir << "'\n";
//        cout << "        (" << n_diff+1 << " mapped potentials  /  md = " << max_square_dst << ")" << endl;
//    }
}
//!============================================================================================================================

//!============================================================================================================================
//!Torsionspotentiale einladen:
void SCORER::load_torsion_potentials() {
    if (!silent_mode) cout << " --> loading torsion potentials" << endl;
    string ptmp = pot_dir + "/torsions.keys";
    ifstream f_in;
    string row,key_name,last_key;
    f_in.open(ptmp.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << ptmp << "'" << endl;
        clean_mem();
        exit(1);
    }
    getline(f_in,row);
    f_in.close();
    istringstream is;
    is.str(row);
    int key_num;
    last_key = "*";
    int n_diff = 0;
    while (!is.eof()) {
        key_name = last_key;
        is >> key_name;
        if (key_name == last_key) break; //!wichtig, um nicht fuer den letzten key 2mal einen eintrag zu machen!!!
        key_num = -1;
        is >> key_num;
        if (key_num == -1) break;
        tors_key_map[key_name] = key_num;
        if (key_num > n_diff) n_diff = key_num;
    }
    string shelp = "/torsion_potentials.bin";
    ptmp = pot_dir + shelp;
    tors_sgrid.load(ptmp.c_str(),true);
//    if (!silent_mode) {
//        cout << "    -> " << tors_key_map.size() << " different torsion potentials loaded from '" << pot_dir << "'\n";
//        cout << "        (" << n_diff+1 << " mapped torsions)" << endl;
//    }
}
//!============================================================================================================================

void SCORER::load_sas_potentials() {
    if (!silent_mode) cout << " --> loading ligand sas potentials" << endl;
    string ptmp = pot_dir + "/sas_potentials.keys";
    ifstream f_in;
    string row,key_name,last_key;
    f_in.open(ptmp.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << ptmp << "'" << endl;
        clean_mem();
        exit(1);
    }
    getline(f_in,row);
    f_in.close();
    istringstream is;
    is.str(row);
    int key_num;
    last_key = "*";
    int n_diff = 0;
    while (!is.eof()) {
        key_name = last_key;
        is >> key_name;
        if (key_name == last_key) break; //!wichtig, um nicht fuer den letzten key 2mal einen eintrag zu machen!!!
        key_num = -1;
        is >> key_num;
        if (key_num == -1) break;
        sas_key_map[key_name] = key_num;
        if (key_num > n_diff) n_diff = key_num;
    }
    string shelp = "/sas_potentials.bin";
    ptmp = pot_dir + shelp;
    sas_grid.load(ptmp.c_str(),true);
//    if (!silent_mode) {
//        cout << "    -> " << sas_key_map.size() << " different ligand sas potentials loaded from '" << pot_dir << "'" << endl;
//    }
}

void SCORER::load_pro_sas_potentials() {
    if (!silent_mode) cout << " --> loading protein sas potentials" << endl;
    string ptmp = pot_dir + "/pro_sas_potentials.keys";
    ifstream f_in;
    string row,key_name,last_key;
    f_in.open(ptmp.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << ptmp << "'" << endl;
        clean_mem();
        exit(1);
    }
    getline(f_in,row);
    f_in.close();
    istringstream is;
    is.str(row);
    int key_num;
    last_key = "*";
    int n_diff = 0;
    while (!is.eof()) {
        key_name = last_key;
        is >> key_name;
        if (key_name == last_key) break; //!wichtig, um nicht fuer den letzten key 2mal einen eintrag zu machen!!!
        key_num = -1;
        is >> key_num;
        if (key_num == -1) break;
        pro_sas_key_map[key_name] = key_num;
        if (key_num > n_diff) n_diff = key_num;
    }
    string shelp = "/pro_sas_potentials.bin";
    ptmp = pot_dir + shelp;
    pro_sas_grid.load(ptmp.c_str(),true);
//    if (!silent_mode) {
//        cout << "    -> " << pro_sas_key_map.size() << " different protein sas potentials loaded from '" << pot_dir << "'" << endl;
//    }
}

//!============================================================================================================================
//!Strukturen einladen:
void SCORER::read_gold_water() {
    string row, hlp;
    vec3d<float> wc;
    //Fuer jeden Liganden einen "Wasserliganden" erzeugen, der als Atoms die O.3 s bekommt, die
    //in den Gold-Solutions als 'ON' markiert sind
    ifstream f_in;
    f_in.open(lig_file.c_str());
    for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
        stl_ptr<LIGAND> wlig = new LIGAND();
        int cid = 1;
        wlig->name = "gold_water_lig";
        while (!f_in.eof()) {
            getline(f_in,row);
            istringstream is;
            is.str(row);
            is >> hlp; is >> hlp;
            if (hlp == "<Gold.Protein.RotatedWaterAtoms>") {
                bool breaker = false;
                while (!f_in.eof()) {
                    getline(f_in,row);
                    if (f_in.eof()) break;
                    if (row[0] == '>') {
                        breaker = true;
                        break;
                    }
                    is.clear(); is.str(row);
                    try {
                        hlp = "###";
                        is >> wc[0] >> wc[1] >> wc[2] >> hlp;
                        if (hlp != "O.3") continue;
                        
                        while (!is.eof()) is >> hlp; //! weil bei GOLD 5 das ON 2 zeichen spaeter kommt
                        //is >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp
                        //>> hlp >> hlp >> hlp >> hlp >> hlp >> hlp >> hlp;
                        
                        if (hlp == "ON") {

                            //cerr << "Wasser marsch: coords = " << wc << endl;

                            stl_ptr<ATOM> watom = new ATOM();
                            watom->intern_type = "O.h2o";
                            watom->coord = wc;
                            
                            watom->name = "O";
                            watom->res_name = "HOH";
                            watom->res_number = 1;
                            watom->intern_id = cid; ++cid;
                            
                            wlig->atoms.push_back(watom);
                            continue;
                        }
                        
                    } catch(...) {
                        continue;
                    }
                }
                if (breaker) break;
            }
        }
        wat_structure->ligands.push_back(wlig);
    }
    f_in.close();
}

//!============================================================================================================================
//!Strukturen einladen:
void SCORER::load_structures() {
    bool kill_exts = true;
    if (calc_rmsd || score_torsions || jiggle) kill_exts = false;
    int uverb = 0;//verbosity;
    //Protein:
    string pro_def = pot_dir + "/protein.def";
    if (!silent_mode) cout << " --> loading protein from " << pro_file << endl;
    pro_structure = new STRUCTURE(verbosity);
    pro_parser = new PARSER(&(*pro_structure),verbosity,false,true,true,false); //Liganden und Wasser mit parsen
    if (pro_type == "mol") {
        if (!(pro_parser->read_mol2(pro_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << pro_file << endl;
            clean_mem();
            exit(1);
        }
        if (!silent_mode) cout << "    -> setting atom types for " << pro_file << endl;
        for (ligands_vec lig=pro_structure->ligands.begin(); lig!=pro_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,uverb,true,true);
        }
    } else {
        if (!(pro_parser->read_pdb(pro_file.c_str()))) {
            cerr << c_message<cERROR>("could not load ") << pro_file << endl;
            clean_mem();
            exit(1);
        }
        
    //    if (pro_structure->ligands.size() > 0) cerr << " got lig " << pro_structure->ligands[0]->name << endl;
        
        if (!pro_structure->protein) {
            cerr << c_message<cWARNING>("found no protein atoms in ") << pro_file << endl;
        }
        
        if (!silent_mode) cout << "    -> setting atom types for " << pro_file << endl;
        
        //pro_structure->pdb2mol2(0,pro_def.c_str(),true,false);
        if (!pro_structure->protein) {
            cerr << c_message<cERROR>("found no protein atoms in ") << pro_file << endl;
            exit(1);
        }
        for (chains_vec cv=pro_structure->protein->chains.begin(); cv!=pro_structure->protein->chains.end(); ++cv) {
            (*cv)->get_atom_typing(1,true,pro_def.c_str(),false,0,true,true);
        }

        //!Jetzt koennen noch Metalle in metals oder ligands stecken:
        for (metals_vec mt=pro_structure->metals.begin(); mt!=pro_structure->metals.end(); ++mt) {
            stl_ptr<LIGAND> mlig = new LIGAND();
            mlig->atoms.push_back((*mt)->atom);
            pro_structure->ligands.push_back(mlig);
        }
        
        pro_structure->metals.clear(); //!DEBUG 190808
        
        for (ligands_vec lig=pro_structure->ligands.begin(); lig!=pro_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,0,true,true);
        }
    }
    //Liganden:
    string lig_def = pot_dir + "/ligand.def";
    if (!silent_mode) cout << " --> loading ligand(s) from " << lig_file << endl;
    lig_structure = new STRUCTURE(verbosity);
    
    //! Wenn minimiert wird wird ja ein neues File rausgeschrieben => wenn da die Wasserstoffe fehlen
    //! sind die IDs anders => in der Visualisierung werden Torsionen falsch definiert:
    lig_parser = new PARSER(&(*lig_structure),verbosity,false,true,true,gold_water);
    
    if (lig_type == "dlg") { //Liganden aus dlg File lesen
        flex_map.clear();
        if (!(lig_parser->read_dlg(flex_map,pro_structure->protein,lig_file.c_str(),flex_res))) {
            cerr << c_message<cERROR>("could not load ") << lig_file << endl;
            clean_mem();
            exit(1);
        }
    } else { //Liganden aus mol2 File lesen
        if (!(lig_parser->read_mol2(lig_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << lig_file << endl;
            clean_mem();
            exit(1);
        }
    }
    if (visualization) {
        max_ligatm = 0;
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end(); ++lig) {
            if (int((*lig)->atoms.size()) > max_ligatm) max_ligatm = (*lig)->atoms.size();
        }
    }
    if (!silent_mode) cout << "    -> setting atom types for " << lig_file << endl;
    for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end(); ++lig) {
        //! 31.07.2010: alle '|' durch '_' ersetzen:
        for (unsigned int i=0; i<(*lig)->name.size(); ++i) if ((*lig)->name[i] == '|') (*lig)->name[i] = '_';
        if (lig_type == "dlg") {
            (*lig)->dlg2mol2(1,lig_def.c_str(),true,1,kill_exts,true);
        } else {
            if (score_torsions || jiggle) (*lig)->get_atom_typing(1,true,lig_def.c_str(),true,uverb,kill_exts,true);
            else if (visualization) (*lig)->get_atom_typing(1,true,lig_def.c_str(),true,uverb,kill_exts,true);
            else (*lig)->get_atom_typing(1,true,lig_def.c_str(),false,uverb,kill_exts,true);
        }
    }
    //Kofaktoren:
    //!Wenn der Kofaktor (gleiches fuer Wasser und Metalle) nicht einfach zum Protein gezaehlt werden soll,
    //!Dann muss er sowohl die Protein-, als auch die Ligand-Atomtypen bekommen. Das entsprechende def_file
    //!wird dann als alt_def uebergeben. Der "alternative" Atomtyp steht dann im Attribut "dict_type"!!!
    string alt_def = "X";
    if (cof_file != "X") {
        if (!silent_mode) cout << " --> loading cofactor(s) from " << cof_file << endl;
        cof_structure = new STRUCTURE(verbosity);
        cof_parser = new PARSER(&(*cof_structure),verbosity,false,true,true,false);
        if (!(cof_parser->read_mol2(cof_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << cof_file << endl;
            clean_mem();
            exit(1);
        }
        if (interaction_mode == 0 || interaction_mode > 4) {
            if (cof_structure->ligands.size() != lig_structure->ligands.size()) {
                cerr << c_message<cERROR>("number of MOLECULE entries in ") << cof_file 
                     << "is different from number of molecules in "
                     << lig_file << ". This is not possible in mode '-I " << interaction_mode << "' !" << endl;
                clean_mem();
                exit(1);
            }
        }
        if (visualization) {
            max_cofatm = 0;
            for (ligands_vec lig=cof_structure->ligands.begin(); lig!=cof_structure->ligands.end(); ++lig) {
                if (int((*lig)->atoms.size()) > max_cofatm) max_cofatm = (*lig)->atoms.size();
            }
        }
        if (interaction_mode == 0 || interaction_mode > 4) alt_def = lig_def;
        if (!silent_mode) cout << "    -> setting atom types for " << cof_file << endl;
        for (ligands_vec lig=cof_structure->ligands.begin(); lig!=cof_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,uverb,true,true,alt_def.c_str());
        }
    }
    //Wasser:
    alt_def = "X";
    if (gold_water) {
        if (!silent_mode) cout << " --> reading information about ON switched water from " << lig_file << endl;
        wat_structure = new STRUCTURE(verbosity);
        read_gold_water();
        if (visualization) {
            max_watatm = 0;
            for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
                if (int((*lig)->atoms.size()) > max_watatm) max_watatm = (*lig)->atoms.size();
            }
        }
        if (interaction_mode == 0 || interaction_mode == 2 ||
            interaction_mode == 4 || interaction_mode == 6) alt_def = lig_def;
        if (!silent_mode) cout << "    -> setting atom types for water molecules" << endl;
        for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,uverb,true,true,alt_def.c_str());
        }
        if (wat_file != "X") {
            cerr << c_message<cWARNING>("You have specified a file with extra water molecules, but you also switched ON ")
                 << "the GOLD modus, where water is read from the ligands file. Your extra file will not be regarded!!!" << endl;
            wat_file = "X";
        }
    } else if (wat_file != "X") {
        if (!silent_mode) cout << " --> loading water(s) from " << wat_file << endl;
        wat_structure = new STRUCTURE(verbosity);
        wat_parser = new PARSER(&(*wat_structure),verbosity,false,true,true,false);
        if (!(wat_parser->read_mol2(wat_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << wat_file << endl;
            clean_mem();
            exit(1);
        }
        if (interaction_mode == 0 || interaction_mode == 2 ||
            interaction_mode == 4 || interaction_mode == 6) {
            if (wat_structure->ligands.size() != lig_structure->ligands.size()) {
                cerr << c_message<cERROR>("number of MOLECULE entries in ") << wat_file 
                     << "is different from number of molecules in "
                     << lig_file << ". This is not possible in mode '-I " << interaction_mode << "' !" << endl;
                clean_mem();
                exit(1);
            }
        }
        if (visualization) {
            max_watatm = 0;
            for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
                if (int((*lig)->atoms.size()) > max_watatm) max_watatm = (*lig)->atoms.size();
            }
        }
        if (interaction_mode == 0 || interaction_mode == 2 ||
            interaction_mode == 4 || interaction_mode == 6) alt_def = lig_def;
        if (!silent_mode) cout << "    -> setting atom types for " << wat_file << endl;
        for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,uverb,true,true,alt_def.c_str());
        }
    }
    //Metalle:
    alt_def = "X";
    if (met_file != "X") {
        if (!silent_mode) cout << " --> loading metal(s) from " << met_file << endl;
        met_structure = new STRUCTURE(verbosity);
        met_parser = new PARSER(&(*met_structure),verbosity,false,true,true,false);
        if (!(met_parser->read_mol2(met_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << met_file << endl;
            clean_mem();
            exit(1);
        }
        if (interaction_mode == 0 || interaction_mode == 3 ||
            interaction_mode == 4 || interaction_mode == 7) {
            if (met_structure->ligands.size() != lig_structure->ligands.size()) {
                cerr << c_message<cERROR>("number of MOLECULE entries in ") << met_file 
                     << "is different from number of molecules in "
                     << lig_file << ". This is not possible in mode '-I " << interaction_mode << "' !" << endl;
                clean_mem();
                exit(1);
            }
        }
        if (visualization) {
            max_metatm = 0;
            for (ligands_vec lig=met_structure->ligands.begin(); lig!=met_structure->ligands.end(); ++lig) {
                if (int((*lig)->atoms.size()) > max_metatm) max_metatm = (*lig)->atoms.size();
            }
        }
        if (interaction_mode == 0 || interaction_mode == 3 ||
            interaction_mode == 4 || interaction_mode == 7) alt_def = lig_def;
        if (!silent_mode) cout << "    -> setting atom types for " << met_file << endl;
        for (ligands_vec lig=met_structure->ligands.begin(); lig!=met_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,uverb,true,true,alt_def.c_str());
        }
    }
    //Referenz fuer rmsd-Berechnung:
    if (ref_file != "X" && calc_rmsd) {
        if (!silent_mode) cout << " --> loading reference for rmsd calculation from " << ref_file << endl;
        ref_structure = new STRUCTURE(verbosity);
        ref_parser = new PARSER(&(*ref_structure),verbosity,false,true,true,false);
        if (!(ref_parser->read_mol2(ref_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << ref_file << endl;
            clean_mem();
            exit(1);
        }
        if (interaction_mode == 0 || interaction_mode == 3 ||
            interaction_mode == 4 || interaction_mode == 7) alt_def = lig_def;
        if (!silent_mode) cout << "    -> setting atom types for " << ref_file << "  for rmsd calculation" << endl;
        for (ligands_vec lig=ref_structure->ligands.begin(); lig!=ref_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,uverb,kill_exts,true,alt_def.c_str());
        }
    }
}
//!============================================================================================================================

//!============================================================================================================================
double SCORER::ENERGY::operator()(VecDoub_I &point) {
    //! 1.) Translation und Rotation:
    vec3d<float> ta;
    ta[0] = point[0];
    ta[1] = point[1];
    ta[2] = point[2];
    matrix<float> rmx(1.,0.,0.,
                      0.,cos(point[3]),-sin(point[3]),
                      0.,sin(point[3]),cos(point[3]));
    
    matrix<float> rmy(cos(point[4]),0,sin(point[4]),
                      0.,1.,0.,
                      -sin(point[4]),0.,cos(point[4]));
    
    matrix<float> rmz(cos(point[5]),-sin(point[5]),0.,
                      sin(point[5]),cos(point[5]),0.,
                      0.,0.,1.);
    
    for (unsigned int i=0; i<lig_structure->ligands[jig_ligno]->atoms.size(); ++i) {
        lig_structure->ligands[jig_ligno]->atoms[i]->coord = ref_coords[i];
        lig_structure->ligands[jig_ligno]->atoms[i]->coord -= centroid;
        lig_structure->ligands[jig_ligno]->atoms[i]->coord *= rmx;
        lig_structure->ligands[jig_ligno]->atoms[i]->coord *= rmy;
        lig_structure->ligands[jig_ligno]->atoms[i]->coord *= rmz;
        lig_structure->ligands[jig_ligno]->atoms[i]->coord += centroid;
        lig_structure->ligands[jig_ligno]->atoms[i]->coord += ta;
    }
    
    //! 2.) Torsionen:
    int pnum = 6;
    for (map<stl_ptr<BOND>,vector<stl_ptr<ATOM> > >::iterator it=r_rot.begin(); it!=r_rot.end(); ++it) {
        vec3d<float> bax(it->first->to->coord);
        bax -= it->first->from->coord;
        bax.norm();
        matrix<float> rm = rotmatrix(bax,point[pnum]); ++pnum;
        for (atoms_vec at=it->second.begin(); at!=it->second.end(); ++at) {
            (*at)->coord -= it->first->to->coord;
            (*at)->coord *= rm;
            (*at)->coord += it->first->to->coord;
        }
    }
    
    //! 3.) TorsionsScore bestimmen:
    double score = 0.;
    double tbuf;
    int sgrid_index[2];
    
    for (map<stl_ptr<BOND>,multimap<string,vector<stl_ptr<ATOM> > > >::iterator it=r_tors.begin(); it!=r_tors.end(); ++it) {
        tbuf = 0.;
        for (multimap<string,vector<stl_ptr<ATOM> > >::iterator ang=it->second.begin(); ang!=it->second.end(); ++ang) {
            float angle = dihedral(ang->second[0]->coord,ang->second[1]->coord,ang->second[2]->coord,ang->second[3]->coord);
            int bin_key = int(0.5+(angle / tors_bin_size));
            if (bin_key > max_tors_bin) bin_key = max_tors_bin;
            sgrid_index[0] = tors_key_map[ang->first];
            sgrid_index[1] = bin_key;
            tbuf += tors_sgrid.value(sgrid_index);
        }
        if (it->second.size()) {
            score += tbuf / it->second.size();
        }
    }
    score *= torsion_weight;
    
    //! 4.) Score bestimmen:
    for (unsigned int i=0; i<lig_structure->ligands[jig_ligno]->atoms.size(); ++i) {
        stl_ptr<ATOM> alig = lig_structure->ligands[jig_ligno]->atoms[i];
        if (alig->sybyl_type == "X") continue;
        for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,alig->coord); apro!=tree_pro->end(); ++apro) {
            float sdst = get_square_distance(alig->coord,(*apro)->coord);
            if (sdst > max_square_dst) continue;
            sdst = sqrt(sdst);
            int bin_key = int(sdst*100.);
            //string ckey = alig->sybyl_type + "_" + (*apro)->sybyl_type;
            string ckey = (*apro)->sybyl_type + "_" + alig->sybyl_type;
            if (key_map.find(ckey) == key_map.end()) continue;
            sgrid_index[0] = key_map[ckey];
            sgrid_index[1] = bin_key;
            score += sgrid.value(sgrid_index);
        }
    }
    
    //! 4.) intramolekulare Clashes:
    float clash_score = 0.;
    for (map<stl_ptr<ATOM>,vector<stl_ptr<ATOM> > >::iterator it=clash_at.begin(); it!=clash_at.end(); ++it) {
        for (atoms_vec at=it->second.begin(); at!=it->second.end(); ++at) {
            float sdst = get_square_distance((*at)->coord,it->first->coord);
            if (sdst > max_square_dst) continue;
            sdst = sqrt(sdst);
            int bin_key = int(sdst*100.);
            string ckey = (*at)->sybyl_type + "_" + it->first->sybyl_type;
            if (key_map.find(ckey) == key_map.end()) continue;
            sgrid_index[0] = key_map[ckey];
            sgrid_index[1] = bin_key;
            if (sgrid.value(sgrid_index) > 0.) clash_score += sgrid.value(sgrid_index);
        }
    }
    score += intra_weight * clash_score;
    
    return score;
}
//!============================================================================================================================

void SCORER::get_r_bonded(stl_ptr<BOND> &to,stl_ptr<ATOM> &curr,tr1::unordered_set<int> &visited) {
    for (atoms_vec at=curr->bonded_atoms.begin(); at!=curr->bonded_atoms.end(); ++at) {
        if (visited.find((*at)->intern_id) != visited.end()) continue;
        visited.insert((*at)->intern_id);
        r_rot[to].push_back(*at);
        get_r_bonded(to,*at,visited);
    }
}

//!============================================================================================================================
//!Und Scoren:
void SCORER::score() {
#ifdef _DSX_DEBUG
    string dhlp1 = pro_file;
    string_fu::remove_ext(dhlp1);
    dhlp1 = string_fu::get_text_after(dhlp1,'/');
    string dhlp2 = lig_file;
    string_fu::remove_ext(dhlp2);
    dhlp2 = string_fu::get_text_after(dhlp2,'/');
    string df1 = "DSX_" + dhlp1 + "_" + dhlp2 + "_debug_scores.txt";
    string df2 = "DSX_" + dhlp1 + "_" + dhlp2 + "_debug_torsions.txt";
    string df3 = "DSX_" + dhlp1 + "_" + dhlp2 + "_debug_matrix.txt";
    cout << " --> running in debug mode: writing pairwise scores to : " << df1 << "\n";
    cout << "                            and torsion scores to      : " << df2 << endl;
    cout << "                            and score matrix to        : " << df3 << endl;
    debug_out.clear();
    debug_out.open(df1.c_str());
    if (score_torsions) {
        debug_tors_out.clear();
        debug_tors_out.open(df2.c_str());
    }
    matrix_out.clear();
    matrix_out.open(df3.c_str());
#endif
    verbosity = 1;
    if (silent_mode) verbosity = 0;
    
    clock_t start,ende; //Variablen zum Zeitstoppen
    double zeit,t1;//
    
    start=clock();
    //! 1.) Potentiale laden:
    load_potentials();
    if (score_torsions || jiggle) load_torsion_potentials();
    if (score_sas) {
        load_sas_potentials();
        load_pro_sas_potentials();
    }
    
    //! 2.) Strukturen einladen:
    load_structures();
    
    //! 3.) Proteinatome in einen Vector 'grid_atoms' packen:
    grid_atoms.clear();
    if (pro_type == "mol") {
        for (ligands_vec lig=pro_structure->ligands.begin(); lig!=pro_structure->ligands.end(); ++lig) {
            for (atoms_vec apro=(*lig)->atoms.begin(); apro!=(*lig)->atoms.end(); ++apro) {
                grid_atoms.push_back(*apro);
                if (as_profile) {
                    if (res_name.find((*apro)->res_number) == res_name.end()) res_name[(*apro)->res_number] = (*apro)->res_name;
                }
            }
        }
    } else {
        if (pro_structure->protein)
        for (chains_vec pch=pro_structure->protein->chains.begin(); pch!=pro_structure->protein->chains.end(); ++pch) {
            for (atoms_vec apro=(*pch)->atoms.begin(); apro!=(*pch)->atoms.end(); ++apro) {
                grid_atoms.push_back(*apro);
                if (as_profile) {
                    if (res_name.find((*apro)->res_number) == res_name.end()) res_name[(*apro)->res_number] = (*apro)->res_name;
                }
            }
        }
    #ifdef _DSX_DEBUG
        debug_out << "--> detected metals in pdb protein (using as part of the protein):" << "\n";
    #endif
        for (ligands_vec lig=pro_structure->ligands.begin(); lig!=pro_structure->ligands.end(); ++lig) {
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                if ((*at)->type == 3) {
                    grid_atoms.push_back(*at); //!Metall
                #ifdef _DSX_DEBUG
                    debug_out << "     " << (*at)->name << " (" << (*at)->sybyl_type << ")" << "\n";
                #endif
                }
            }
        }
    }
    
    if (jiggle) {
        if (cof_file != "X") for (ligands_vec cof=cof_structure->ligands.begin(); cof!=cof_structure->ligands.end(); ++cof) {
            for (atoms_vec acof=(*cof)->atoms.begin(); acof!=(*cof)->atoms.end(); ++acof) {
                grid_atoms.push_back(*acof);
            }
        }
        if (wat_file != "X") for (ligands_vec cof=wat_structure->ligands.begin(); cof!=wat_structure->ligands.end(); ++cof) {
            for (atoms_vec acof=(*cof)->atoms.begin(); acof!=(*cof)->atoms.end(); ++acof) {
                grid_atoms.push_back(*acof);
            }
        }
        if (met_file != "X") for (ligands_vec cof=met_structure->ligands.begin(); cof!=met_structure->ligands.end(); ++cof) {
            for (atoms_vec acof=(*cof)->atoms.begin(); acof!=(*cof)->atoms.end(); ++acof) {
                grid_atoms.push_back(*acof);
            }
        }
    }
    
    int n_ligands = lig_structure->ligands.size();

    tree_pro = new OCTREE<ATOM>(grid_atoms,max_dst);

    if (jiggle) {
        if (!silent_mode) cout << " --> minimizing ligands" << endl;
        if (n_ligands < 5) spacing = 3.;//4.;
        else if (n_ligands < 20) spacing = 2.;//3.;
        else if (n_ligands < 100) spacing = 1.;//2.;
        else spacing = 0.5;
        hsp = sqrt(3. * spacing * spacing);
        hsp = (hsp/2.) + 0.2;
        max_lim = max_dst + hsp;
        max_lim *= max_lim;

        ENERGY energy;
        
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
            ref_coords.clear();
            vec3d<float> np(0.,0.,0.);
            for (atoms_vec at=lig_structure->ligands[jig_ligno]->atoms.begin();
                           at!=lig_structure->ligands[jig_ligno]->atoms.end(); ++at) {
                ref_coords.push_back((*at)->coord);
                np += (*at)->coord;
            }
            np /= float(ref_coords.size());
            centroid = np;
            
            int frb = (*lig)->get_freerot_bonds();
            r_rot.clear();
            for (bonds_vec bt=(*lig)->bonds.begin(); bt!=(*lig)->bonds.end();++bt) {
                if (!(*bt)->free_rot) continue;
                r_rot[*bt] = vector<stl_ptr<ATOM> >();
                tr1::unordered_set<int> visited; visited.insert((*bt)->to->intern_id); visited.insert((*bt)->from->intern_id);
                get_r_bonded(*bt,(*bt)->to,visited);
            }
            
            r_tors.clear();
            for (map<stl_ptr<BOND>,vector<stl_ptr<ATOM> > >::iterator bt=r_rot.begin(); bt!=r_rot.end(); ++bt) {
                if (bt->first->from->ext->n_heavy_bonded < 2 || bt->first->to->ext->n_heavy_bonded < 2) continue;
                r_tors[bt->first] = multimap<string,vector<stl_ptr<ATOM> > >();
                for (atoms_vec bfa=bt->first->from->bonded_atoms.begin(); bfa!=bt->first->from->bonded_atoms.end(); ++bfa) {
                    if ((*bfa)->intern_id == bt->first->to->intern_id) continue;
                    for (atoms_vec bta=bt->first->to->bonded_atoms.begin(); bta!=bt->first->to->bonded_atoms.end(); ++bta) {
                        if ((*bta)->intern_id == bt->first->from->intern_id) continue;
                        
                        string adder = "1"; //! weder 'from' noch 'to' ist Teil eines Ringes
                        if (bt->first->from->ext->is_ring) {
                            if (bt->first->to->ext->is_ring) {
                                adder = "3"; //! from und to gehoeren zu 2 verschiedenen Ringen
                                stl_ptr<RING> fring,tring;
                                for (rings_vec rit=lig_structure->ligands[jig_ligno]->rings.begin();
                                rit!=lig_structure->ligands[jig_ligno]->rings.end(); ++rit) {
                                    bool ffrom = false;
                                    bool fto = false;
                                    for (atoms_vec at=(*rit)->ring_atoms.begin();
                                            at!=(*rit)->ring_atoms.end(); ++at) {
                                        if (*at == bt->first->from) {
                                            ffrom = true;
                                            fring = *rit;
                                        }
                                        if (*at == bt->first->to) {
                                            fto = true;
                                            tring = *rit;
                                        }
                                    }
                                    if (ffrom && fto) { // beide im gleichen Ring
                                        adder = "4"; //! from und to gehoeren zum selben Ring
                                        break;
                                    }
                                }
                                if (adder == "3") {
                                    int n_cat = 0;
                                    for (atoms_vec at=fring->ring_atoms.begin();
                                            at!=fring->ring_atoms.end(); ++at) {
                                        for (atoms_vec bt=tring->ring_atoms.begin();
                                                bt!=tring->ring_atoms.end(); ++bt) {
                                            if (*at == *bt) {++n_cat; break;}
                                        }
                                        if (n_cat > 1) break;
                                    }
                                    if (n_cat > 1) adder = "4";
                                }
                            } else adder = "2"; //! nur from ODER to ist Teil eines Ringes
                        } else if (bt->first->to->ext->is_ring) adder = "2"; //! nur from ODER to ist Teil eines Ringes
                        
                        string ckey = (*bfa)->intern_type + "_" + bt->first->from->intern_type + "_" +
                                bt->first->to->intern_type + "_" + (*bta)->intern_type + "_" + adder;
                        
                        if (tors_key_map.find(ckey) == tors_key_map.end()) {
                            if (!try_to_replace_with_reduced((*bfa)->intern_type,bt->first->from->intern_type,
                                                             bt->first->to->intern_type,(*bta)->intern_type,
                                                             adder,ckey)) continue;
                        }

                        multimap<string,vector<stl_ptr<ATOM> > >::iterator my;
                        my = r_tors[bt->first].insert(pair<string,vector<stl_ptr<ATOM> > >(ckey,vector<stl_ptr<ATOM> >()));
                        my->second.push_back(*bfa);
                        my->second.push_back(bt->first->from);
                        my->second.push_back(bt->first->to);
                        my->second.push_back(*bta);
                    }
                }
            }
            
            clash_at.clear();
            if ((*lig)->SP_map == 0) (*lig)->calc_SP_map();
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                if ((*at)->sybyl_type == "X") continue;
                atoms_vec posi = at; ++posi;
                if (posi == (*lig)->atoms.end()) break;
                clash_at[*at] = vector<stl_ptr<ATOM> >();
                for (atoms_vec bt=posi; bt!=(*lig)->atoms.end(); ++bt) {
                    if ((*bt)->sybyl_type == "X") continue;
                    if ((*bt)->intern_id == (*at)->intern_id) continue;
                    if ((*lig)->SP_map[(*at)->intern_id][(*bt)->intern_id] > 3) clash_at[*at].push_back(*bt);
                }
            }
            
            freedom = 6 + frb;
            double* ap = new double[freedom];
            for (int i=0; i<freedom; ++i) {
                ap[i] = 0.;
            }
            const VecDoub_I point(freedom,ap);
            
            Powell<ENERGY> powell(energy,10.);

            try {
                VecDoub pmin = powell.minimize(point);
            } catch(...) {
                cout << "warning: exceeded max number of iterations in line search" << endl;
            }
            
            delete[] ap;
            ++jig_ligno;

            if (jiggle_first_only) {
                if (jig_ligno == 5) break;
            }
        }
        
        string jname = lig_file;
        string_fu::replace_ext(jname,".mol2","_jiggled.mol2");
        lig_parser->write_mol2(jname.c_str());
        lig_file = jname;
        
        //! Jetzt grid_atoms neu erstellen:
        grid_atoms.clear();
        if (pro_type == "mol") {
            for (ligands_vec lig=pro_structure->ligands.begin(); lig!=pro_structure->ligands.end(); ++lig) {
                for (atoms_vec apro=(*lig)->atoms.begin(); apro!=(*lig)->atoms.end(); ++apro) {
                    grid_atoms.push_back(*apro);
                    if (as_profile) {
                        if (res_name.find((*apro)->res_number) == res_name.end()) 
                            res_name[(*apro)->res_number] = (*apro)->res_name;
                    }
                }
            }
        } else {
            if (pro_structure->protein)
            for (chains_vec pch=pro_structure->protein->chains.begin(); pch!=pro_structure->protein->chains.end(); ++pch) {
                for (atoms_vec apro=(*pch)->atoms.begin(); apro!=(*pch)->atoms.end(); ++apro) {
                    grid_atoms.push_back(*apro);
                    if (as_profile) {
                        if (res_name.find((*apro)->res_number) == res_name.end()) 
                            res_name[(*apro)->res_number] = (*apro)->res_name;
                    }
                }
            }
            for (ligands_vec lig=pro_structure->ligands.begin(); lig!=pro_structure->ligands.end(); ++lig) {
                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                    if ((*at)->type == 3) grid_atoms.push_back(*at); //!Metall
                }
            }
        }
    }
    
    if (!silent_mode && jiggle) {
        ende=clock();
        zeit=(ende-start);
        t1 = zeit/CLOCKS_PER_SEC;
        cout << "Minimizing the ligands  : " << t1 << " s" << endl;
        start=clock();
    }
    
    //! 4.) Protein-Grid vorbereiten:
    if (n_ligands < 3) spacing = 10.;
    else if (n_ligands < 6) spacing = 8.;
    else if (n_ligands < 12) spacing = 6.;
    else if (n_ligands < 24) spacing = 4.;
    else if (n_ligands < 100) spacing = 2.;
    //else if (n_ligands < 1000) spacing = 2.;
    else spacing = 1.;
    /*
    if (n_ligands < 3) spacing = 12.;
    else if (n_ligands < 10) spacing = 8.;
    else if (n_ligands < 30) spacing = 5.;
    else if (n_ligands < 100) spacing = 4.;
    else if (n_ligands < 1000) spacing = 3.;
    else if (n_ligands < 10000) spacing = 2.;
    else spacing = 1.;
    */    
    hsp = sqrt(3. * spacing * spacing);
    hsp = (hsp/2.) + 0.2;
    max_lim = max_dst + hsp;
    max_lim *= max_lim;
    
    int sizebuf[2];
    int sizebuf2[2];
    int sgrid_index[2];
    if (visualization) {
        sizebuf[0] = lig_structure->ligands.size();
        sizebuf[1] = grid_atoms.size();
        int vis_indi = 0;
        for (atoms_vec at=grid_atoms.begin(); at!=grid_atoms.end(); ++at) {
            (*at)->id = vis_indi;
            ++vis_indi;
        }
        pro_vis.resize(sizebuf,0.);
        sizebuf[1] = max_ligatm;
        lig_vis.resize(sizebuf,0.);
        if (cof_file != "X") {
            sizebuf[1] = max_cofatm;
            cof_vis.resize(sizebuf,0.);
        }
        if (wat_file != "X" || gold_water) {
            sizebuf[1] = max_watatm;
            wat_vis.resize(sizebuf,0.);
        }
        if (met_file != "X") {
            sizebuf[1] = max_metatm;
            met_vis.resize(sizebuf,0.);
        }
        for (unsigned int i=0; i<lig_structure->ligands.size(); ++i) {
            good_con[i] = vector<TRIPLE>();
            bad_con[i] = vector<TRIPLE>();
            if (score_torsions) bad_tors[i] = vector<int>();
        }
    }

    if (!silent_mode) {
        ende=clock();//
        zeit=(ende-start);//
        t1 = zeit/CLOCKS_PER_SEC;//
        cout << "Preparing the structures  : " << t1 << " s" << endl;
        start=clock();
    }
    
    if (!silent_mode) start=clock();//
    //! 5.) Jetzt das eigentliche Scoring:
    int bin_key;
    int ligno = 0;
    int n_contacts;
    int n_tors = 0;
    int per_tors;
    float score;
    float tors_score;
    float intra_score;
    float sas_score;
    float sdst;
    float angle;
    float tbuf;
    string ckey;

    if (score_sas) {
        vector<stl_ptr<ATOM> > sas_protein_atoms;
        for (atoms_vec alit=grid_atoms.begin(); alit!=grid_atoms.end(); ++alit) {
            if (pro_sas_key_map.find((*alit)->intern_type) == sas_key_map.end()) continue;
            sas_protein_atoms.push_back(*alit);
        }
        sas_tree_pro = new OCTREE<ATOM>(sas_protein_atoms,MAX_RAD2,MAX_RAD2);
    }

    if (tree_pro) {
        delete tree_pro;
        tree_pro = new OCTREE<ATOM>(grid_atoms,max_dst);
    }

    if (!silent_mode) cout << " --> calculating scores for " << lig_structure->ligands.size() << " structures" << endl;
    
    #ifdef _DSX_DEBUG
        debug_out << "--> atom-atom-scores fuer " << lig_structure->ligands.size() << " Liganden:" << "\n";
        debug_out << "    (Bei Interaktionsmodus Null werden nur die Interaktionen zum Liganden," << "\n";
        debug_out << "     NICHT zum Protein beruecksichtigt)" << "\n";
    #endif
    
    //Ueber alle Liganden:
    for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
        score = 0.;
        tors_score = 0.;
        intra_score = 0.;
        sas_score = 0.;
        n_contacts = 0;
        sizebuf[0] = ligno;
        sizebuf2[0] = ligno;
        
        if (as_profile) {
            res_score[ligno] = map<int,float>();
            for (map<int,string>::iterator rt=res_name.begin(); rt!=res_name.end(); ++rt) {
                res_score[ligno][rt->first] = 0.;
            }
        }
        
        if (score_torsions) {
            #ifdef _DSX_DEBUG
            debug_tors_out << "\n==============================================================================================\n";
            debug_tors_out << "==> " << (*lig)->name << ":" << "\n";
            debug_tors_out << "==============================================================================================\n";
            debug_tors_out << "   Number   |             Torsion             |     Score     | Free rotatable | No Potential\n";
            debug_tors_out << "------------|---------------------------------|---------------|----------------|--------------\n";
            map<int,int> tors_num;
            int curr_tors_num = 0;
            #endif
            n_tors = 0;
            for (bonds_vec bt=(*lig)->bonds.begin(); bt!=(*lig)->bonds.end();++bt) {
                if ((*bt)->from->ext->n_heavy_bonded < 2 || (*bt)->to->ext->n_heavy_bonded < 2) continue;
                if (!((*bt)->from->element == "C" || (*bt)->from->element == "N" || (*bt)->from->element == "O" ||
				      (*bt)->from->element == "S") || (*bt)->from->element == "P") continue;
				if (!((*bt)->to->element == "C" || (*bt)->to->element == "N" || (*bt)->to->element == "O" ||
				      (*bt)->to->element == "S") || (*bt)->to->element == "P") continue;
                tbuf = 0.;
                per_tors = 0;
                for (atoms_vec bfa=(*bt)->from->bonded_atoms.begin(); bfa!=(*bt)->from->bonded_atoms.end(); ++bfa) {
                    if ((*bfa)->intern_id == (*bt)->to->intern_id) continue;
                    //if ((*bfa)->intern_type[0] == 'H') continue;
                    for (atoms_vec bta=(*bt)->to->bonded_atoms.begin(); bta!=(*bt)->to->bonded_atoms.end(); ++bta) {
                        if ((*bta)->intern_id == (*bt)->from->intern_id) continue;
                        //if ((*bta)->intern_type[0] == 'H') continue;
                        string adder = "1"; //! weder 'from' noch 'to' ist Teil eines Ringes
                        if ((*bt)->from->ext->is_ring) {
                            if ((*bt)->to->ext->is_ring) {
                                adder = "3"; //! from und to gehoeren zu 2 verschiedenen Ringen
                                stl_ptr<RING> fring,tring;
                                for (rings_vec rit=(*lig)->rings.begin(); rit!=(*lig)->rings.end(); ++rit) {
                                    bool ffrom = false;
                                    bool fto = false;
                                    for (atoms_vec at=(*rit)->ring_atoms.begin();
                                                   at!=(*rit)->ring_atoms.end(); ++at) {
                                        if (*at == (*bt)->from) {
                                            ffrom = true;
                                            fring = *rit;
                                        }
                                        if (*at == (*bt)->to) {
                                            fto = true;
                                            tring = *rit;
                                        }
                                    }
                                    if (ffrom && fto) { // beide im gleichen Ring
                                        adder = "4"; //! from und to gehoeren zum selben Ring
                                        break;
                                    }
                                }
                                if (adder == "3") {
                                    int n_cat = 0;
                                    for (atoms_vec at=fring->ring_atoms.begin();
                                                   at!=fring->ring_atoms.end(); ++at) {
                                        for (atoms_vec bt=tring->ring_atoms.begin();
                                                           bt!=tring->ring_atoms.end(); ++bt) {
                                            if (*at == *bt) {++n_cat; break;}
                                        }
                                        if (n_cat > 1) break;
                                    }
                                    if (n_cat > 1) adder = "4";
                                }
                            } else adder = "2"; //! nur from ODER to ist Teil eines Ringes
                        } else if ((*bt)->to->ext->is_ring) adder = "2"; //! nur from ODER to ist Teil eines Ringes
                        
                        ckey = (*bfa)->intern_type + "_" + (*bt)->from->intern_type + "_" +
                               (*bt)->to->intern_type + "_" + (*bta)->intern_type + "_" + adder;
                        
                        #ifdef _DSX_DEBUG
                        string debug_key;
                        if ((*bta)->intern_type < (*bfa)->intern_type) {
                            debug_key = (*bta)->intern_type + "_" + (*bt)->to->intern_type + "_" +
                                            (*bt)->from->intern_type + "_" + (*bfa)->intern_type + "_" + adder;
                        } else {
                            debug_key = (*bfa)->intern_type + "_" + (*bt)->from->intern_type + "_" +
                                            (*bt)->to->intern_type + "_" + (*bta)->intern_type + "_" + adder;
                        }
                        int tors_ident;
                        if ((*bt)->to->intern_id < (*bt)->from->intern_id) {
                            tors_ident = 100000 * (*bt)->to->intern_id + (*bt)->from->intern_id;
                        } else {
                            tors_ident = 100000 * (*bt)->from->intern_id + (*bt)->to->intern_id;
                        }
                        if (tors_num.find(tors_ident) == tors_num.end()) {
                            tors_num[tors_ident] = curr_tors_num;
                            ++curr_tors_num;
                        }
                        if (tors_key_map.find(ckey) == tors_key_map.end()) {
                            debug_tors_out << " ";
                            debug_tors_out.width(10); debug_tors_out << tors_num[tors_ident] << " | ";
                            debug_tors_out.width(31); debug_tors_out << debug_key << " | ";
                            debug_tors_out.width(13); debug_tors_out << 0.0 << " | ";
                            debug_tors_out.width(14); 
                            if ((*bt)->free_rot) debug_tors_out << "yes" << " | ";
                            else debug_tors_out << "no" << " | ";
                            debug_tors_out << "yes\n";
                        }
                        #endif
                        
                        if (tors_key_map.find(ckey) == tors_key_map.end()) {
                            if (!try_to_replace_with_reduced((*bfa)->intern_type,(*bt)->from->intern_type,
                                                             (*bt)->to->intern_type,(*bta)->intern_type,
                                                             adder,ckey)) {
                                                             continue;
                            }
                        }
                        
                        angle = dihedral((*bfa)->coord,(*bt)->from->coord,(*bt)->to->coord,(*bta)->coord);
                        
                        bin_key = int(0.5+(angle / tors_bin_size));
                        if (bin_key > max_tors_bin) bin_key = max_tors_bin;
                        
                        sgrid_index[0] = tors_key_map[ckey];
                        sgrid_index[1] = bin_key;
                        
                        #ifdef _DSX_DEBUG
                        debug_tors_out << " ";
                        debug_tors_out.width(10); debug_tors_out << tors_num[tors_ident] << " | ";
                        debug_tors_out.width(31); debug_tors_out << debug_key << " | ";
                        float debug_val = tors_sgrid.value(sgrid_index);
                        debug_tors_out.width(13);debug_tors_out << debug_val << " | ";
                        debug_tors_out.width(14);
                        if ((*bt)->free_rot) debug_tors_out << "yes" << " | ";
                        else debug_tors_out << "no" << " | ";
                        debug_tors_out << "no\n";
                        #endif
                        
                        ///if (!(*bt)->free_rot) {
                        ///    if (tors_sgrid.value(sgrid_index) < 0.) continue;
                        ///}
                        
                        tbuf += tors_sgrid.value(sgrid_index);
                        ++per_tors;
                        
                        if (visualization) {
                            if (tors_sgrid.value(sgrid_index) > vis_tors) {
                                bad_tors[ligno].push_back((*bfa)->id);
                                bad_tors[ligno].push_back((*bt)->from->id);
                                bad_tors[ligno].push_back((*bt)->to->id);
                                bad_tors[ligno].push_back((*bta)->id);
                            }
                        }
                    }
                }
                if (per_tors > 0) {
                    tors_score += tbuf / per_tors;
                    ++n_tors;
                }
            }
        }

        if (score_intra) {
            if (!jiggle) {
                if ((*lig)->SP_map == 0) (*lig)->calc_SP_map();
            }
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                atoms_vec posi = at; ++posi;
                for (atoms_vec bt=posi; bt!=(*lig)->atoms.end(); ++bt) {
                    if ((*lig)->SP_map[(*at)->intern_id][(*bt)->intern_id] > 3) {
                        float sdst = get_square_distance((*at)->coord,(*bt)->coord);
                        if (sdst > max_square_dst) continue;
                        sdst = sqrt(sdst);
                        int bin_key = int(sdst*100.);
                        string ckey = (*at)->sybyl_type + "_" + (*bt)->sybyl_type;
                        if (key_map.find(ckey) == key_map.end()) continue;
                        sgrid_index[0] = key_map[ckey];
                        sgrid_index[1] = bin_key;
                        if (sgrid.value(sgrid_index) > 0.) intra_score += sgrid.value(sgrid_index);
                    }
                }
            }
        }

        if (score_sas) {
            vector<SAS_POINT*> lig_sas_points;

            tr1::unordered_map<ATOM*,uint64_t> t_frei; //! ratio

            tr1::unordered_map<ATOM*,uint64_t> t_komp;
            tr1::unordered_set<ATOM*> pro_atoms_to_check;
            vector<stl_ptr<ATOM> > ligand_atoms;
            for (atoms_vec alit=(*lig)->atoms.begin(); alit!=(*lig)->atoms.end(); ++alit) {
                if (sas_key_map.find((*alit)->intern_type) == sas_key_map.end()) continue;
                ligand_atoms.push_back(*alit);
            }
            OCTREE<ATOM> sas_tree_lig(ligand_atoms,LMAX_RAD2,LMAX_RAD2);
            for (atoms_vec alit=ligand_atoms.begin(); alit!=ligand_atoms.end(); ++alit) {
                t_komp[alit->get_pointer()] = 0;
                t_frei[alit->get_pointer()] = 0;
                for (vector<vec3d<float> >::iterator st=sphere_grid_coords::surf_coords[(*alit)->element].begin();
                                                     st!=sphere_grid_coords::surf_coords[(*alit)->element].end(); ++st) {
                    SAS_POINT* cp = new SAS_POINT((*st) + (*alit)->coord);
                    cp->closest_lig = alit->get_pointer();
                    cp->lig_dst = sphere_grid_coords::sas_square_dst[(*alit)->element];
                    for (vector<ATOM*>::iterator alit2=sas_tree_lig.begin(sas_lig_lig_spacing,cp->coord); alit2!=sas_tree_lig.end(); ++alit2) {
                        if ((*alit)->intern_id == (*alit2)->intern_id) continue;
                        float sdst = get_square_distance((*alit2)->coord,cp->coord);
                        if (sdst < sphere_grid_coords::sas_square_dst[(*alit2)->element]) {
                            //! Entweder nicht frei oder anderem lig_atom zugehoeriger Oberflaechenteil
                            delete cp;
                            cp = 0;
                            break;
                        } else if (sdst < cp->lig_dst) {
                            cp->lig_dst = sdst;
                            cp->closest_lig = *alit2;
                        }
                    }
                    if (cp) lig_sas_points.push_back(cp);
                }
            }
            for (vector<SAS_POINT*>::iterator st=lig_sas_points.begin(); st!=lig_sas_points.end(); ++st) {
                for (vector<ATOM*>::iterator apit=sas_tree_pro->begin(sas_lig_pro_spacing,(*st)->coord); apit!=sas_tree_pro->end(); ++apit) {
                    float sdst = get_square_distance((*apit)->coord,(*st)->coord);
                    if (sdst < sphere_grid_coords::sas_square_dst[(*apit)->element]) {
                        (*st)->is_opposite = true;
                        if (sdst < (*st)->pro_dst) {
                            (*st)->pro_dst = sdst;
                            (*st)->closest_pro = *apit; //wird nur gebraucht, wenn is_opposite
                        }
                    }
                }
                if ((*st)->is_opposite) {
                    if ((*st)->closest_pro->element == "O" || (*st)->closest_pro->element == "N") {
                        if ((*st)->closest_lig->element == "O" || (*st)->closest_lig->element == "N") {
                            (*st)->is_opposite = false;
                        }
                    }
                    pro_atoms_to_check.insert((*st)->closest_pro);
                }
                t_frei[(*st)->closest_lig] += 1;
                if ((*st)->is_opposite == false) {
                    t_komp[(*st)->closest_lig] += 1;
                }
            }
            for (tr1::unordered_map<ATOM*,uint64_t>::iterator ct=t_komp.begin(); ct!=t_komp.end(); ++ct) {
                sgrid_index[0] = sas_key_map[ct->first->intern_type];

                #if defined (MODE_RATIO)
                if (t_frei[ct->first] == 0) t_frei[ct->first] = 1;
                sgrid_index[1] = 100 * ct->second / t_frei[ct->first];
                #elif defined (MODE_DELTA)
                sgrid_index[1] = t_frei[ct->first] - ct->second;
                #endif

                sas_score += sas_grid.value(sgrid_index);
            }
            for (vector<SAS_POINT*>::iterator st=lig_sas_points.begin(); st!=lig_sas_points.end(); ++st) delete *st;



            vector<SAS_POINT*> pro_sas_points;
            tr1::unordered_map<ATOM*,uint64_t> pro_t_komp;
            tr1::unordered_map<ATOM*,uint64_t> pro_t_frei;
            for (tr1::unordered_set<ATOM*>::iterator alit=pro_atoms_to_check.begin(); alit!=pro_atoms_to_check.end(); ++alit) {
                pro_t_komp[*alit] = 0;
                pro_t_frei[*alit] = 0;
                for (vector<vec3d<float> >::iterator st=sphere_grid_coords::surf_coords[(*alit)->element].begin();
                                                     st!=sphere_grid_coords::surf_coords[(*alit)->element].end(); ++st) {
                    SAS_POINT* cp = new SAS_POINT((*st) + (*alit)->coord);
                    cp->closest_lig = *alit;
                    cp->lig_dst = sphere_grid_coords::sas_square_dst[(*alit)->element];
                    for (vector<ATOM*>::iterator alit2=sas_tree_pro->begin(sas_lig_lig_spacing,cp->coord); alit2!=sas_tree_pro->end(); ++alit2) {
                        if ((*alit)->intern_id == (*alit2)->intern_id) continue;
                        float sdst = get_square_distance((*alit2)->coord,cp->coord);
                        if (sdst < sphere_grid_coords::sas_square_dst[(*alit2)->element]) {
                            //! Entweder nicht frei oder anderem lig_atom zugehoeriger Oberflaechenteil
                            delete cp;
                            cp = 0;
                            break;
                        } else if (sdst < cp->lig_dst) {
                            cp->lig_dst = sdst;
                            cp->closest_lig = *alit2;
                        }
                    }
                    if (cp) pro_sas_points.push_back(cp);
                }
            }
            for (vector<SAS_POINT*>::iterator st=pro_sas_points.begin(); st!=pro_sas_points.end(); ++st) {
                for (vector<ATOM*>::iterator apit=sas_tree_lig.begin(sas_lig_lig_spacing,(*st)->coord); apit!=sas_tree_lig.end(); ++apit) {
                    float sdst = get_square_distance((*apit)->coord,(*st)->coord);
                    if (sdst < sphere_grid_coords::sas_square_dst[(*apit)->element]) {
                        (*st)->is_opposite = true;
                        if (sdst < (*st)->pro_dst) {
                            (*st)->pro_dst = sdst;
                            (*st)->closest_pro = *apit; //wird nur gebraucht, wenn is_opposite
                        }
                    }
                }
                if ((*st)->is_opposite) {
                    if ((*st)->closest_pro->element == "O" || (*st)->closest_pro->element == "N") {
                        if ((*st)->closest_lig->element == "O" || (*st)->closest_lig->element == "N") {
                            (*st)->is_opposite = false;
                        }
                    }
                }
                pro_t_frei[(*st)->closest_lig] += 1;
                if ((*st)->is_opposite == false) {
                    pro_t_komp[(*st)->closest_lig] += 1;
                }
            }
            for (tr1::unordered_map<ATOM*,uint64_t>::iterator ct=pro_t_komp.begin(); ct!=pro_t_komp.end(); ++ct) {
                sgrid_index[0] = pro_sas_key_map[ct->first->intern_type];


                #if defined (MODE_RATIO)
                if (pro_t_frei[ct->first] == 0) pro_t_frei[ct->first] = 1;
                sgrid_index[1] = 100 * ct->second / pro_t_frei[ct->first];
                #elif defined (MODE_DELTA)
                sgrid_index[1] = pro_t_frei[ct->first] - ct->second;
                #endif
                
                sas_score += pro_sas_grid.value(sgrid_index);
            }
            for (vector<SAS_POINT*>::iterator st=pro_sas_points.begin(); st!=pro_sas_points.end(); ++st) delete *st;
        }

        if (coval_check) {
            //Alle Protein-Ligand-Kombinationen durchgehen und Abstaende bestimmen.
            //Bei Ueberlagerung und gleichem Element wird von einer kovalenten Bindung ausgegangen.
            //Die Sybyltypen beider Atome werden auf "X" gesetzt, ebenso wie die Typen der jeweils
            //direkt gebundenen Atome. => Die Atome werden im Scoring dann nicht beruecksichtigt.
            //Bevor zum naechstem Komplex uebergegangen wird muessen die Typen wieder zurueckgesetzt werden!!!
            //=> Urspruenglichen Sybyltyp in dem Attribut dict_type zwischenlagern.
            for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                (*alig)->dict_type = (*alig)->sybyl_type;
                for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*alig)->coord); apro!=tree_pro->end(); ++apro) {
                    (*apro)->dict_type = (*apro)->sybyl_type;
                }
            }
            for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                sizebuf2[1] = alig - (*lig)->atoms.begin();
                for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*alig)->coord); apro!=tree_pro->end(); ++apro) {
                    if ((*apro)->sybyl_type[0] != (*alig)->sybyl_type[0]) continue;
                    if (flex_res) {
                        if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                            sdst = get_square_distance((*alig)->coord,(*apro)->coord);
                        } else sdst = get_square_distance((*alig)->coord,flex_map[ligno][(*apro)->intern_id]);
                    } else sdst = get_square_distance((*alig)->coord,(*apro)->coord);
                    if (sdst < 0.035) {
                        (*alig)->sybyl_type = "X";
                        (*apro)->sybyl_type = "X";
                        for (atoms_vec bat=(*alig)->bonded_atoms.begin();
                                       bat!=(*alig)->bonded_atoms.end(); ++bat) {
                            (*bat)->sybyl_type = "X";
                        }
                        if (pro_type == "mol") {
                            for (atoms_vec bat=(*apro)->bonded_atoms.begin();
                                       bat!=(*apro)->bonded_atoms.end(); ++bat) {
                                (*bat)->sybyl_type = "X";
                            }
                        } else { //es gibt keine bonded_atoms
                            for (vector<ATOM*>::iterator api=tree_pro->begin(spacing,(*alig)->coord); api!=tree_pro->end(); ++api) {
                                vec3d<float> hav,hbv;
                                if (flex_res) {
                                    if (flex_map[ligno].find((*apro)->intern_id) ==
                                        flex_map[ligno].end()) {
                                        hav = (*apro)->coord;
                                    } else hav = flex_map[ligno][(*apro)->intern_id];
                                    if (flex_map[ligno].find((*api)->intern_id) ==
                                        flex_map[ligno].end()) {
                                        hbv = (*api)->coord;
                                    } else hbv = flex_map[ligno][(*api)->intern_id];
                                } else {
                                    hav = (*apro)->coord;
                                    hbv = (*api)->coord;
                                }
                                sdst = get_square_distance(hav,hbv);
                                if (sdst < 4.) (*api)->sybyl_type = "X";
                            }
                        }
                    }
                }
            }
        }
        
    #ifdef _DSX_DEBUG
        debug_out << "\n" << "=========================================================================================" << "\n";
        debug_out << "==> " << (*lig)->name << ":" << "\n";
        debug_out << "\n" << "=================================================================================" << "\n";
        debug_out << "***                       Protein-Ligand Wechselwirkungen                     ***" << "\n";
        debug_out << "=================================================================================" << "\n";
        debug_out << "--> Einzelbeitraege zum Score:" << "\n";
        debug_out << "     Ligand atom | Receptor atom        | L Type | R Type | distance | PPI score " << "\n";
        debug_out << "    -----------------------------------------------------------------------------" << "\n";
        
        //! Matrix erstellen, die alle Pro-Lig-Interaktionen enthaelt
        //!     rows    = Ligandatome
        //!     columns = Proteinatome
        matrix_out << "\n# new ligand: " << (*lig)->name << "\n";
        set<int> col_ids;
        vector<int> row_ids;
        for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
            if ((*alig)->sybyl_type == "X") continue;
            row_ids.push_back((*alig)->intern_id);
            sizebuf2[1] = alig - (*lig)->atoms.begin();
            for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*alig)->coord); apro!=tree_pro->end(); ++apro) {
                col_ids.insert((*apro)->intern_id);
            }
        }
        matrix_out << "n_columns = " << col_ids.size() << "\n";
        matrix_out << "n_rows = " << row_ids.size() << "\n";
        matrix_out << "columns: ";
        tr1::unordered_map<int,int> id2num;
        int cnum = 0;
        for (set<int>::iterator cit=col_ids.begin(); cit!=col_ids.end(); ++cit) {
            matrix_out << *cit << " ";
            id2num[*cit] = cnum; ++cnum;
        }
        matrix_out << "\n";
        matrix_out << "rows: ";
        for (vector<int>::iterator cit=row_ids.begin(); cit!=row_ids.end(); ++cit) {
            matrix_out << *cit << " ";
        }
        matrix_out << "\n";
    #endif
        
        //Ueber alle Ligandatome:
        if (score_pairs) {
        for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
            
            if ((*alig)->sybyl_type == "X") continue;

            #ifdef _DSX_DEBUG
            matrix_out << "\n";
            int last_pro = *(col_ids.begin());
            #endif

            sizebuf2[1] = alig - (*lig)->atoms.begin();
            for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*alig)->coord); apro!=tree_pro->end(); ++apro) {
                if (flex_res) {
                    if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                        sdst = get_square_distance((*alig)->coord,(*apro)->coord);
                    } else sdst = get_square_distance((*alig)->coord,flex_map[ligno][(*apro)->intern_id]);
                } else sdst = get_square_distance((*alig)->coord,(*apro)->coord);

                if (sdst > max_square_dst) continue;

                sdst = sqrt(sdst);
                bin_key = int(sdst*100.);

                //ckey = (*alig)->sybyl_type + "_" + (*apro)->sybyl_type;
                ckey = (*apro)->sybyl_type + "_" + (*alig)->sybyl_type;

                if (key_map.find(ckey) == key_map.end()) {
                    //cerr << "NOKEY " << ckey << endl;
                    continue;
                }

                ++n_contacts;

                sgrid_index[0] = key_map[ckey];
                sgrid_index[1] = bin_key;

                tbuf = sgrid.value(sgrid_index);
                score += tbuf;

                if (as_profile) {
                    res_score[ligno][(*apro)->res_number] += tbuf;
                }
            #ifdef _DSX_DEBUG
                debug_out << "     ";
                debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";

                debug_out.width(7); debug_out << right << (*apro)->name; debug_out << " (";
                debug_out.width(3); debug_out << (*apro)->res_name; debug_out << " " << (*apro)->chain_id;
                debug_out << " "; debug_out.width(4); debug_out << (*apro)->res_number; debug_out << ") | ";

                debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";

                debug_out.width(6); debug_out << (*apro)->sybyl_type; debug_out << " | ";

                debug_out.width(8); debug_out << sdst; debug_out << " | ";

                debug_out << tbuf << "\n";

                int cdiff = id2num[(*apro)->intern_id] - id2num[last_pro];
                for (int ci=1; ci<cdiff; ++ci) matrix_out << "/ ";
                matrix_out << tbuf << " ";
                last_pro = (*apro)->intern_id;
            #endif
                if (visualization) {
                    sizebuf[1] = (*apro)->id;
                    pro_vis.value(sizebuf) += tbuf;
                    lig_vis.value(sizebuf2) += tbuf;
                    if (sdst < specific_inter_dst) {
                        if (tbuf < vis_lower) {
                            if (flex_res) {
                                if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*apro)->coord,tbuf));
                                } else good_con[ligno].push_back(TRIPLE((*alig)->coord,
                                       flex_map[ligno][(*apro)->intern_id],tbuf));
                            } else good_con[ligno].push_back(TRIPLE((*alig)->coord,(*apro)->coord,tbuf));
                        } else if (tbuf > vis_upper) {
                            if (flex_res) {
                                if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*apro)->coord,tbuf));
                                } else bad_con[ligno].push_back(TRIPLE((*alig)->coord,
                                       flex_map[ligno][(*apro)->intern_id],tbuf));
                            } else bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*apro)->coord,tbuf));
                        }
                    }
                }
            }
        }
    
        //Nochmal fuer Interaktionen mit dem Kofaktor:
        if (cof_file != "X") {
            
        #ifdef _DSX_DEBUG
            debug_out << "\n" << "=================================================================================" << "\n";
            debug_out << "***                       Ligand-Cofactor Wechselwirkungen                    ***" << "\n";
            debug_out << "=================================================================================" << "\n";
            debug_out << "--> Einzelbeitraege zum Score:" << "\n";
            debug_out << "     Ligand atom | Cofactor atom | L Type | R Type | distance | PPI score " << "\n";
            debug_out << "    ----------------------------------------------------------------------" << "\n";
        #endif
        
            if (interaction_mode > 0 && interaction_mode < 5) { //Kofaktor als Teil des Proteins:
                for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                sizebuf2[1] = alig - (*lig)->atoms.begin();
                //Ueber alle Kofaktoren:
                for (ligands_vec cof=cof_structure->ligands.begin(); cof!=cof_structure->ligands.end(); ++cof) {
                    //Ueber alle Atome des aktuellen Kofaktors
                    for (atoms_vec acof=(*cof)->atoms.begin(); acof!=(*cof)->atoms.end(); ++acof) {
                        sdst = get_square_distance((*alig)->coord,(*acof)->coord);
                        
                        if (sdst > max_square_dst) continue;
                        
                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);
                
                        //ckey = (*alig)->sybyl_type + "_" + (*acof)->sybyl_type;
                        ckey = (*acof)->sybyl_type + "_" + (*alig)->sybyl_type;

                        if (key_map.find(ckey) == key_map.end()) continue;
                        
                        ++n_contacts;
                
                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;
                        
                    #ifdef _DSX_DEBUG
                        debug_out << "     ";
                        debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";
                        
                        debug_out.width(5); debug_out << right << (*acof)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*acof)->id; debug_out << ")   | ";
                        
                        debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";
                        
                        debug_out.width(6); debug_out << (*acof)->sybyl_type; debug_out << " | ";
                        
                        debug_out.width(8); debug_out << sdst; debug_out << " | ";
                        
                        debug_out << tbuf << "\n";
                    #endif
                
                        if (visualization) {
                            sizebuf[1] = acof - (*cof)->atoms.begin();
                            cof_vis.value(sizebuf) += tbuf;
                            lig_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*acof)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*acof)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
                }
            } else { //Kofaktor individuell behandeln:
                //Hier muss jetzt die Interaktion 'Ligand Nr.X mit Kofaktor Nr.X', sowie die
                //Interaktion 'Protein mit Kofaktor Nr.X' gescored werden.
                for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                    sizebuf2[1] = alig - (*lig)->atoms.begin();
                    for (atoms_vec acof=cof_structure->ligands[ligno]->atoms.begin();
                                   acof!=cof_structure->ligands[ligno]->atoms.end(); ++acof) { //Ligand<-->Kofaktor
                        sdst = get_square_distance((*alig)->coord,(*acof)->coord);

                        if (sdst > max_square_dst) continue;

                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);

                        //ckey = (*alig)->sybyl_type + "_" + (*acof)->sybyl_type;
                        ckey = (*acof)->sybyl_type + "_" + (*alig)->sybyl_type;

                        if (key_map.find(ckey) == key_map.end()) continue;

                        ++n_contacts;

                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;

                    #ifdef _DSX_DEBUG
                        debug_out << "     ";
                        debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";

                        debug_out.width(5); debug_out << right << (*acof)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*acof)->id; debug_out << ") | ";

                        debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";

                        debug_out.width(6); debug_out << (*acof)->sybyl_type; debug_out << " | ";

                        debug_out.width(8); debug_out << sdst; debug_out << " | ";

                        debug_out << tbuf << "\n";
                    #endif

                        if (visualization) {
                            sizebuf[1] = acof - cof_structure->ligands[ligno]->atoms.begin();
                            cof_vis.value(sizebuf) += tbuf;
                            lig_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*acof)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*acof)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
        
                for (atoms_vec acof=cof_structure->ligands[ligno]->atoms.begin(); 
                               acof!=cof_structure->ligands[ligno]->atoms.end(); ++acof) { //Kofaktor<-->Protein
                    //Uber alle relevanten Proteinatome:
                    for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*acof)->coord); apro!=tree_pro->end(); ++apro) {
                        if (flex_res) {
                            if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                                sdst = get_square_distance((*acof)->coord,(*apro)->coord);
                            } else sdst =
                              get_square_distance((*acof)->coord,flex_map[ligno][(*apro)->intern_id]);
                        } else sdst = get_square_distance((*acof)->coord,(*apro)->coord);

                        if (sdst > max_square_dst) continue;

                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);

                        //ckey = (*acof)->dict_type + "_" + (*apro)->sybyl_type;
                        ckey = (*apro)->sybyl_type + "_" + (*acof)->dict_type;

                        if (key_map.find(ckey) == key_map.end()) continue;

                        ++n_contacts;

                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;

                        if (as_profile) {
                            res_score[ligno][(*apro)->res_number] += tbuf;
                        }

                        if (visualization) {
                            sizebuf[1] = (*apro)->id;
                            pro_vis.value(sizebuf) += tbuf;
                            sizebuf2[1] = acof - cof_structure->ligands[ligno]->atoms.begin();
                            cof_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    if (flex_res) {
                                        if (flex_map[ligno].find((*apro)->intern_id) ==
                                                                 flex_map[ligno].end()) {
                                            good_con[ligno].push_back(TRIPLE((*acof)->coord,
                                                                      (*apro)->coord,tbuf));
                                        } else good_con[ligno].push_back(TRIPLE((*acof)->coord,
                                               flex_map[ligno][(*apro)->intern_id],tbuf));
                                    } else good_con[ligno].push_back(TRIPLE((*acof)->coord,
                                                                     (*apro)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    if (flex_res) {
                                        if (flex_map[ligno].find((*apro)->intern_id) ==
                                                                 flex_map[ligno].end()) {
                                            bad_con[ligno].push_back(TRIPLE((*acof)->coord,
                                                                     (*apro)->coord,tbuf));
                                        } else bad_con[ligno].push_back(TRIPLE((*acof)->coord,
                                               flex_map[ligno][(*apro)->intern_id],tbuf));
                                    } else bad_con[ligno].push_back(TRIPLE((*acof)->coord,
                                                                    (*apro)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
            }
        }
    
        //Nochmal fuer Interaktionen mit Wasser:
        if (wat_file != "X" || gold_water) {
            
        #ifdef _DSX_DEBUG
            debug_out << "\n" << "=================================================================================" << "\n";
            debug_out << "***                       Ligand-Wasser  Wechselwirkungen                     ***" << "\n";
            debug_out << "=================================================================================" << "\n";
            debug_out << "--> Einzelbeitraege zum Score:" << "\n";
            debug_out << "     Ligand atom | Water atom  | L Type | R Type | distance | PPI score " << "\n";
            debug_out << "    --------------------------------------------------------------------" << "\n";
        #endif
        
            if (interaction_mode == 1 || interaction_mode == 3 ||
                interaction_mode == 5 || interaction_mode == 7) { //Wasser als Teil des Proteins:
                for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                sizebuf2[1] = alig - (*lig)->atoms.begin();
                //Ueber alle Kofaktoren:
                int watno = -1;
                for (ligands_vec wat=wat_structure->ligands.begin(); wat!=wat_structure->ligands.end(); ++wat) {
                    ++watno;
                    if (gold_water && ligno != watno) continue;
                    //Ueber alle Atome der aktuellen Wasseratome
                    for (atoms_vec awat=(*wat)->atoms.begin(); awat!=(*wat)->atoms.end(); ++awat) {
                        sdst = get_square_distance((*alig)->coord,(*awat)->coord);
                        
                        if (sdst > max_square_dst) continue;
                        
                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);
                
                        //ckey = (*alig)->sybyl_type + "_" + (*awat)->sybyl_type;
                        ckey = (*awat)->sybyl_type + "_" + (*alig)->sybyl_type;
                        
                        if (key_map.find(ckey) == key_map.end()) continue;
                        
                        ++n_contacts;
                
                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;
                            
                    #ifdef _DSX_DEBUG
                        debug_out << "     ";
                        debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";
                        
                        debug_out.width(5); debug_out << right << (*awat)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*awat)->id; debug_out << ") | ";
                        
                        debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";
                        
                        debug_out.width(6); debug_out << (*awat)->sybyl_type; debug_out << " | ";
                        
                        debug_out.width(8); debug_out << sdst; debug_out << " | ";
                        
                        debug_out << tbuf << "\n";
                    #endif
                    
                        if (visualization) {
                            sizebuf[1] = awat - (*wat)->atoms.begin();
                            wat_vis.value(sizebuf) += tbuf;
                            lig_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*awat)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*awat)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
                }
            } else { //Wasser individuell behandeln:
                //Hier muss jetzt die Interaktion 'Ligand Nr.X mit Wassercluster Nr.X', sowie die
                //Interaktion 'Protein mit Wassercluster Nr.X' gescored werden.
                for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                sizebuf2[1] = alig - (*lig)->atoms.begin();
                    for (atoms_vec awat=wat_structure->ligands[ligno]->atoms.begin();
                                   awat!=wat_structure->ligands[ligno]->atoms.end(); ++awat) { //Ligand<-->Wasser
                        sdst = get_square_distance((*alig)->coord,(*awat)->coord);

                        if (sdst > max_square_dst) continue;

                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);

                        ckey = (*alig)->sybyl_type + "_" + (*awat)->sybyl_type;

                        if (key_map.find(ckey) == key_map.end()) continue;

                        ++n_contacts;

                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;

                    #ifdef _DSX_DEBUG
                        debug_out << "     ";
                        debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";

                        debug_out.width(5); debug_out << right << (*awat)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*awat)->id; debug_out << ") | ";

                        debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";

                        debug_out.width(6); debug_out << (*awat)->sybyl_type; debug_out << " | ";

                        debug_out.width(8); debug_out << sdst; debug_out << " | ";

                        debug_out << tbuf << "\n";
                    #endif

                        if (visualization) {
                            sizebuf[1] = awat - wat_structure->ligands[ligno]->atoms.begin();
                            wat_vis.value(sizebuf) += tbuf;
                            lig_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*awat)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*awat)->coord,tbuf));
                                }
                            }
                        }
                    }
                }

                for (atoms_vec awat=wat_structure->ligands[ligno]->atoms.begin(); 
                               awat!=wat_structure->ligands[ligno]->atoms.end(); ++awat) { //Wasser<-->Protein
                    for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*awat)->coord); apro!=tree_pro->end(); ++apro) {
                        if (flex_res) {
                            if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                                sdst = get_square_distance((*awat)->coord,(*apro)->coord);
                            } else sdst =
                                   get_square_distance((*awat)->coord,flex_map[ligno][(*apro)->intern_id]);
                        } else sdst = get_square_distance((*awat)->coord,(*apro)->coord);


                        if (sdst > max_square_dst) continue;

                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);

                        //ckey = (*awat)->dict_type + "_" + (*apro)->sybyl_type;
                        ckey = (*apro)->sybyl_type + "_" + (*awat)->dict_type;

                        if (key_map.find(ckey) == key_map.end()) continue;

                        ++n_contacts;

                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;

                        if (as_profile) {
                            res_score[ligno][(*apro)->res_number] += tbuf;
                        }

                        if (visualization) {
                            sizebuf[1] = (*apro)->id;
                            pro_vis.value(sizebuf) += tbuf;
                            sizebuf2[1] = awat -wat_structure->ligands[ligno]->atoms.begin() ;
                            wat_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    if (flex_res) {
                                        if (flex_map[ligno].find((*apro)->intern_id) ==
                                                                 flex_map[ligno].end()) {
                                            good_con[ligno].push_back(TRIPLE((*awat)->coord,
                                            (*apro)->coord,tbuf));
                                        } else good_con[ligno].push_back(TRIPLE((*awat)->coord,
                                               flex_map[ligno][(*apro)->intern_id],tbuf));
                                    } else good_con[ligno].push_back(TRIPLE((*awat)->coord,
                                           (*apro)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    if (flex_res) {
                                        if (flex_map[ligno].find((*apro)->intern_id) ==
                                                                 flex_map[ligno].end()) {
                                            bad_con[ligno].push_back(TRIPLE((*awat)->coord,
                                            (*apro)->coord,tbuf));
                                        } else bad_con[ligno].push_back(TRIPLE((*awat)->coord,
                                               flex_map[ligno][(*apro)->intern_id],tbuf));
                                    } else bad_con[ligno].push_back(TRIPLE((*awat)->coord,
                                           (*apro)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
            }
        }
    
        //Nochmal fuer Interaktionen mit Metallen:
        if (met_file != "X") {
            
        #ifdef _DSX_DEBUG
            debug_out << "\n" << "=================================================================================" << "\n";
            debug_out << "***                       Ligand-Metall  Wechselwirkungen                     ***" << "\n";
            debug_out << "=================================================================================" << "\n";
            debug_out << "--> Einzelbeitraege zum Score:" << "\n";
            debug_out << "     Ligand atom | Metal atom  | L Type | R Type | distance | PPI score " << "\n";
            debug_out << "    --------------------------------------------------------------------" << "\n";
        #endif
        
            if (interaction_mode == 1 || interaction_mode == 2 ||
                interaction_mode == 5 || interaction_mode == 6) { //Metall als Teil des Proteins:
                for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                sizebuf2[1] = alig - (*lig)->atoms.begin();
                //Ueber alle Metalle:
                for (ligands_vec met=met_structure->ligands.begin(); met!=met_structure->ligands.end(); ++met) {
                    //Ueber alle Metallatome
                    for (atoms_vec amet=(*met)->atoms.begin(); amet!=(*met)->atoms.end(); ++amet) {
                        sdst = get_square_distance((*alig)->coord,(*amet)->coord);
                        
                        if (sdst > max_square_dst) continue;
                        
                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);
                
                        ckey = (*alig)->sybyl_type + "_" + (*amet)->sybyl_type;
                        
                        if (key_map.find(ckey) == key_map.end()) continue;
                        
                        ++n_contacts;
                
                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;
                        
                    #ifdef _DSX_DEBUG
                        debug_out << "     ";
                        debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";
                        
                        debug_out.width(5); debug_out << right << (*amet)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*amet)->id; debug_out << ") | ";
                        
                        debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";
                        
                        debug_out.width(6); debug_out << (*amet)->sybyl_type; debug_out << " | ";
                        
                        debug_out.width(8); debug_out << sdst; debug_out << " | ";
                        
                        debug_out << tbuf << "\n";
                    #endif
                        
                        if (visualization) {
                            sizebuf[1] = amet - (*met)->atoms.begin();
                            met_vis.value(sizebuf) += tbuf;
                            lig_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*amet)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*amet)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
                }
            } else { //Metalle individuell behandeln:
                //Hier muss jetzt die Interaktion 'Ligand Nr.X mit Metall Nr.X', sowie die
                //Interaktion 'Protein mit Metall Nr.X' gescored werden.
                for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                    sizebuf2[1] = alig - (*lig)->atoms.begin();
                    for (atoms_vec amet=met_structure->ligands[ligno]->atoms.begin();
                                   amet!=met_structure->ligands[ligno]->atoms.end(); ++amet) { //Ligand<-->Metall
                        sdst = get_square_distance((*alig)->coord,(*amet)->coord);

                        if (sdst > max_square_dst) continue;

                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);

                        ckey = (*alig)->sybyl_type + "_" + (*amet)->sybyl_type;

                        if (key_map.find(ckey) == key_map.end()) continue;

                        ++n_contacts;

                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;

                    #ifdef _DSX_DEBUG
                        debug_out << "     ";
                        debug_out.width(5); debug_out << right << (*alig)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*alig)->id; debug_out << ") | ";

                        debug_out.width(5); debug_out << right << (*amet)->name; debug_out << " (";
                        debug_out.width(3); debug_out << (*amet)->id; debug_out << ") | ";

                        debug_out.width(6); debug_out << (*alig)->sybyl_type; debug_out << " | ";

                        debug_out.width(6); debug_out << (*amet)->sybyl_type; debug_out << " | ";

                        debug_out.width(8); debug_out << sdst; debug_out << " | ";

                        debug_out << tbuf << "\n";
                    #endif

                        if (visualization) {
                            sizebuf[1] = amet - met_structure->ligands[ligno]->atoms.begin();
                            met_vis.value(sizebuf) += tbuf;
                            lig_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    good_con[ligno].push_back(TRIPLE((*alig)->coord,(*amet)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    bad_con[ligno].push_back(TRIPLE((*alig)->coord,(*amet)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
                
                for (atoms_vec amet=met_structure->ligands[ligno]->atoms.begin(); 
                               amet!=met_structure->ligands[ligno]->atoms.end(); ++amet) { //Metall<-->Protein
                    for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*amet)->coord); apro!=tree_pro->end(); ++apro) {
                        if (flex_res) {
                            if (flex_map[ligno].find((*apro)->intern_id) == flex_map[ligno].end()) {
                                sdst = get_square_distance((*amet)->coord,(*apro)->coord);
                            } else sdst =
                                   get_square_distance((*amet)->coord,flex_map[ligno][(*apro)->intern_id]);
                        } else sdst = get_square_distance((*amet)->coord,(*apro)->coord);

                        if (sdst > max_square_dst) continue;

                        sdst = sqrt(sdst);
                        bin_key = int(sdst*100.);

                        //ckey = (*amet)->dict_type + "_" + (*apro)->sybyl_type;
                        ckey = (*apro)->sybyl_type + "_" + (*amet)->dict_type;

                        if (key_map.find(ckey) == key_map.end()) continue;

                        ++n_contacts;

                        tbuf = sgrid[key_map[ckey]][bin_key];
                        score += tbuf;

                        if (as_profile) {
                            res_score[ligno][(*apro)->res_number] += tbuf;
                        }

                        if (visualization) {
                            sizebuf[1] = (*apro)->id;
                            pro_vis.value(sizebuf) += tbuf;
                            sizebuf2[1] = amet - met_structure->ligands[ligno]->atoms.begin();
                            met_vis.value(sizebuf2) += tbuf;
                            if (sdst < specific_inter_dst) {
                                if (tbuf < vis_lower) {
                                    if (flex_res) {
                                        if (flex_map[ligno].find((*apro)->intern_id) ==
                                                                 flex_map[ligno].end()) {
                                            good_con[ligno].push_back(TRIPLE((*amet)->coord,
                                            (*apro)->coord,tbuf));
                                        } else good_con[ligno].push_back(TRIPLE((*amet)->coord,
                                               flex_map[ligno][(*apro)->intern_id],tbuf));
                                    } else good_con[ligno].push_back(TRIPLE((*amet)->coord,
                                           (*apro)->coord,tbuf));
                                } else if (tbuf > vis_upper) {
                                    if (flex_res) {
                                        if (flex_map[ligno].find((*apro)->intern_id) ==
                                                                 flex_map[ligno].end()) {
                                            bad_con[ligno].push_back(TRIPLE((*amet)->coord,
                                            (*apro)->coord,tbuf));
                                        } else bad_con[ligno].push_back(TRIPLE((*amet)->coord,
                                               flex_map[ligno][(*apro)->intern_id],tbuf));
                                    } else bad_con[ligno].push_back(TRIPLE((*amet)->coord,
                                           (*apro)->coord,tbuf));
                                }
                            }
                        }
                    }
                }
            }
        }
        }
        score *= pair_weight;
        if (score_torsions) {
            tors_score_map[ligno] = torsion_weight * tors_score;
            if (n_tors > 0) per_tors_score[ligno] = torsion_weight * tors_score / n_tors;
            score += torsion_weight * tors_score;
        }
        if (score_intra) {
            score += intra_weight * intra_score;
        }
        if (score_sas) {
            sas_score_map[ligno] = sas_weight * sas_score;
            score = score + sas_weight * sas_score;
        }
        score_map[ligno] = score; //! = Score + Torsion_Score
        if ((*lig)->atoms.size() > 0) per_atom_score[ligno] = score / (*lig)->atoms.size();
        else per_atom_score[ligno] = score;
        if (n_contacts == 0) ++n_contacts;
        per_contact_score[ligno] = score / n_contacts;
        ligno++;
        
        if (coval_check) {
            //Die urspruenglichen Sybyltypen wieder herstellen:
            for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                (*alig)->sybyl_type = (*alig)->dict_type;
                for (vector<ATOM*>::iterator apro=tree_pro->begin(spacing,(*alig)->coord); apro!=tree_pro->end(); ++apro) {
                    (*apro)->sybyl_type = (*apro)->dict_type;
                }
            }
        }
    }
    
    if (!silent_mode) {
        ende=clock();//
        zeit=(ende-start);//
        t1 = zeit/CLOCKS_PER_SEC;//
        cout << "Calculation of the scores";
        if (visualization) cout << " and visualization elements";
        cout << "  : " << t1 << " s" << endl;
    }
    
    if ((!silent_mode) && calc_rmsd) start=clock();//
    //! 6.) Gegebenenfalls rmsd-Werte berechnen:
    if (calc_rmsd) {
        if (!silent_mode) cout << " --> calculating rmsd values" << endl;
        ligno = 0;
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
            if (ref_file != "X") {
                rmsd_map[ligno] = ref_structure->ligands[0]->get_rmsd(*lig,false,false);
                if (rmsd_map[ligno] < 0.) rmsd_map[ligno] = ref_structure->ligands[0]->get_bk_rmsd(*lig,false,false);
            } else {
                rmsd_map[ligno] = lig_structure->ligands[0]->get_rmsd(*lig,false,false);
                if (rmsd_map[ligno] < 0.) rmsd_map[ligno] = lig_structure->ligands[0]->get_bk_rmsd(*lig,false,false);
            }
            ligno++;
        }
    }
    
    if ((!silent_mode) && calc_rmsd) {
        ende=clock();//
        zeit=(ende-start);//
        t1 = zeit/CLOCKS_PER_SEC;//
        cout << "Calculation of the rmsd values  : " << t1 << " s" << endl;
    }
    if (!silent_mode) {
        if (tors_key_map.size() > 1) cout << "Final size of tors_key_map: " << tors_key_map.size() << endl;
    }
#ifdef _DSX_DEBUG
    debug_out.close();
    if (score_torsions) debug_tors_out.close();
#endif
}
//!============================================================================================================================

//!============================================================================================================================
//!Ergebnisse in Files schreiben:
void SCORER::write_results(int n_visi,bool write_txt) {
    string hlp1 = pro_file;
    string_fu::remove_ext(hlp1);
    hlp1 = string_fu::get_text_after(hlp1,'/');
    string hlp2 = lig_file;
    string_fu::remove_ext(hlp2);
    hlp2 = string_fu::get_text_after(hlp2,'/');
    string filename = "DSX_" + hlp1 + "_" + hlp2 + ".txt";
    
    if (!silent_mode) cout << " --> writing results to  " << filename << endl;
    
    ofstream f_out;
    
    if (write_txt) {
    f_out.open(filename.c_str());

    f_out << "scoring with 'DSX' version " << version_string << "\n" << "\n";
    f_out << "protein file                  : " << pro_file << "\n";
    f_out << "ligands file                  : " << lig_file << "\n";
    
    f_out << "cofactor file                 : "; if (cof_file == "X") f_out << "none" << "\n"; else f_out << cof_file << "\n";
    f_out << "waters file                   : "; if (wat_file == "X") f_out << "none" << "\n"; else f_out << wat_file << "\n";
    f_out << "metals file                   : "; if (met_file == "X") f_out << "none" << "\n"; else f_out << met_file << "\n";
    
    f_out << "reference file                : "; if (ref_file == "X") f_out << "none" << "\n"; else f_out << ref_file << "\n";
    
    f_out << "used potentials               : " << pot_dir << "\n";
    f_out << "used interaction mode         : " << interaction_mode << "\n";
    f_out << "used sort mode                : " << sort_mode << "\n";
    f_out << "score atom atom pairs         : ";
    if (score_pairs) f_out << "yes" << "\n";
    else f_out << "no" << "\n";
    f_out << "score torsion angles          : ";
    if (score_torsions) f_out << "yes" << "\n";
    else f_out << "no" << "\n";
    f_out << "score intramolecular clashes  : ";
    if (score_intra) f_out << "yes" << "\n";
    else f_out << "no" << "\n";
    f_out << "score solvent access. surf.   : ";
    if (score_sas) f_out << "yes" << "\n";
    else f_out << "no" << "\n";
    f_out << "minimize ligands              : ";
    if (jiggle) f_out << "yes" << "\n";
    else f_out << "no" << "\n";
    f_out << "flexible residues             : ";
    if (flex_res) f_out << "yes" << "\n";
    else f_out << "no" << "\n";
    f_out << "covalent bond check           : ";
    if (coval_check) f_out << "yes" << "\n" << "\n";
    else f_out << "no" << "\n" << "\n";
    f_out << "important notes:" << "\n";
    f_out << "   - The field 'score' is the total score including possible torsion, sas and intramolecular contributions." << "\n";
    f_out << "   - The 'PCS'(per_contact_score) is the score divided by the number of atom-atom-interactions having any" << "\n";
    f_out << "      contribution to the total score (number of contacts within 6A). (Do not confuse with a per atom score.)" << "\n";
    f_out << "   - The 'tors_score' is the sum of scores for each bond. A single bond B--C can have more than one" << "\n";
    f_out << "      torsions (A1--B--C--D1, A2--B--C--D1, ...). The score for a single bond is the mean of its" << "\n";
    f_out << "      possible torsions." << "\n";
    f_out << "   - The 'sas_score' is the solvent accessable surface score for solvation/desolvation contributions." << "\n" << "\n";
    f_out << "@RESULTS" << "\n" << "\n";
    f_out << "  number  |              name              |  rmsd  |   score   |   rank   |    PCS    | tors_score | sas_score \n";
    f_out << "----------|--------------------------------|--------|-----------|----------|-----------|------------|-----------\n";
    }
    
    multimap<float,int> rssorted; //score : number
    for (map<int,float>::iterator it=score_map.begin(); it!=score_map.end(); ++it) {
        rssorted.insert(make_pair(it->second,it->first));
    }
    map<int,int> rank_map; //lig-number:rang
    int rank = 1;
    for (multimap<float,int>::iterator it=rssorted.begin(); it!=rssorted.end(); ++it) {
        rank_map[it->second] = rank;
        ++rank;
    }
    rssorted.clear();
    
    if (write_txt) {
    f_out.setf(ios::fixed,ios::floatfield);
    f_out.precision(3);
    if (sort_mode == 0) { //Originalreihenfolge
        string lname,tester;
        for (map<int,float>::iterator it=score_map.begin(); it!=score_map.end(); ++it) {
            f_out << " "; f_out.width(8); f_out << left << it->first << " | ";
            tester = lig_structure->ligands[it->first]->name;
            if (tester.size() > 30) lname.assign(tester,0,29);
            else if (tester[tester.size()-1] < 33) lname.assign(tester,0,tester.size()-1);
            else lname = tester;
            f_out.width(30); f_out << lname << " | ";
            f_out.width(6); if (calc_rmsd) f_out << rmsd_map[it->first] << " | ";
            else f_out << " none  | ";
            f_out.width(9); f_out << it->second << " | ";
            f_out.width(8); f_out << rank_map[it->first] << " | ";
            f_out.width(9); f_out << per_contact_score[it->first] << " | ";
            f_out.width(10); f_out << tors_score_map[it->first] << " | ";
            f_out.width(9); f_out << sas_score_map[it->first] << "\n";
        }
        f_out << endl;
    } else {
        multimap<float,int> ssorted; //score : number   der gleiche Score kann mehrfach vorkommen => Multimap
        if (sort_mode == 1) for (map<int,float>::iterator it=score_map.begin(); it!=score_map.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        } else if (sort_mode == 2) for (map<int,float>::iterator it=per_atom_score.begin(); it!=per_atom_score.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        } else if (sort_mode == 3) for (map<int,float>::iterator it=per_contact_score.begin(); it!=per_contact_score.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        } else if (sort_mode == 4) for (map<int,float>::iterator it=rmsd_map.begin(); it!=rmsd_map.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        } else if (sort_mode == 5) for (map<int,float>::iterator it=tors_score_map.begin(); it!=tors_score_map.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        } else if (sort_mode == 6) for (map<int,float>::iterator it=per_tors_score.begin(); it!=per_tors_score.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        }
        string lname,tester;
        for (multimap<float,int>::iterator it=ssorted.begin(); it!=ssorted.end(); ++it) {
            f_out << " "; f_out.width(8); f_out << left << it->second << " | ";
            tester = lig_structure->ligands[it->second]->name;
            if (tester.size() > 30) lname.assign(tester,0,29);
            else if (tester[tester.size()-1] < 33) lname.assign(tester,0,tester.size()-1); //!um Steuerzeichen abzuschneiden die sonst
            else lname = tester;                                                           //!die Formatierung zerhauen wuerden
            f_out.width(30); f_out << lname << " | ";
            f_out.width(6); if (calc_rmsd) f_out << rmsd_map[it->second] << " | ";
            else f_out << " none  | ";
            f_out.width(9); f_out << score_map[it->second] << " | ";
            f_out.width(8); f_out << rank_map[it->second] << " | ";
            f_out.width(9); f_out << per_contact_score[it->second] << " | ";
            f_out.width(10); f_out << tors_score_map[it->second] << " | ";
            f_out.width(9); f_out << sas_score_map[it->second] << "\n";
        }
        f_out << endl;
    }
    f_out.close();
    f_out.clear();
    }
    
    if (as_profile) {
        string pname = filename;
        string_fu::replace_ext(pname,".txt","_profiles.txt");
        if (!silent_mode) cout << " --> writing aminoacid profiles to " << pname << endl;
        ofstream p_out;
        p_out.open(pname.c_str());
        
        for (map<int,map<int,float> >::iterator it=res_score.begin(); it!=res_score.end(); ++it) {
            p_out << "--> ligand " << it->first << "\n";
            for (map<int,float>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                if (fabs(jt->second) < 100.) continue;
                p_out << res_name[jt->first] << " " << jt->first << " " << jt->second << "\n";
            }
        }
        
        p_out.close();
    }
    
    //!Und jetzt die Pymol-Visualisierung:
    if (visualization && n_visi == 99999) {
        string wname = "X";
        //!neu 29.07.2010:  Tobi will die GOLD-water auch als mol2 rausgeschrieben haben: ==========================
        if (gold_water) {
            vector<stl_ptr<LIGAND> > ltw;
            for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                    (*at)->sybyl_type = "O.3";
                }
                ltw.push_back(*lig);
            }
            wname = lig_file;
            string_fu::replace_ext(wname,".mol2","_gold_water.mol2");
            if (!silent_mode) cout << " --> writing " << wname << " for visualization" << endl;
            lig_parser->write_mol2(ltw,wname.c_str());
        }
        //!=========================================================================================================
        
        
        string fname = filename;
        string newmolname = lig_file;
        if (flex_res) {
            //Die Cavities als pdb mit verschiedenen MODEL Eintraegen rausschreiben:
            //Die Aenderung der Koordinaten ist an dieser Stelle nur zum Bestimmen der
            //Tasche. In der Funktion write_pdb_cav_all() muessen dann nochmals alle Koordinaten
            //geaendert werden.
            string_fu::replace_ext(fname,".txt","_cavs.pdb");
            if (!silent_mode) cout << " --> flexible residues => writing pdb file with all cavities to  " << fname << endl;
            
            int ligno = 0;
            
            map<int,vec3d<float> > orig_pos;
            for (atoms_vec at=grid_atoms.begin(); at!=grid_atoms.end(); ++at) {
                orig_pos[(*at)->intern_id] = (*at)->coord;
            }
            
            for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                for (atoms_vec at=grid_atoms.begin(); at!=grid_atoms.end(); ++at) {
                    (*at)->coord = orig_pos[(*at)->intern_id];
                    if (flex_map[ligno].find((*at)->intern_id) == flex_map[ligno].end()) continue;
                    (*at)->coord = flex_map[ligno][(*at)->intern_id];
                }
                
                pro_structure->get_lig_cavity(flex_cav_rad,&(**lig));
                ++ligno;
            }
            for (atoms_vec at=grid_atoms.begin(); at!=grid_atoms.end(); ++at) {
                (*at)->coord = orig_pos[(*at)->intern_id];
            }
            pro_parser->write_pdb_cav_all(flex_map,fname.c_str());
            
            if (lig_type == "dlg") {
                string_fu::replace_ext(newmolname,"dlg","mol2");
                if (!silent_mode) cout << " --> dlg input => writing mol2 file for visualization to  " << newmolname << endl;
                lig_parser->write_mol2(newmolname.c_str());
            }
        }
        
        string_fu::replace_ext(filename,"txt","py");
        f_out.open(filename.c_str());
        
        if (!silent_mode) cout << " --> writing pymol visualization file to  " << filename << endl;
        
        f_out << "# This visualization file was automatically created by 'DSX'" << "\n";
        f_out << "# Start pymol and use 'run this_files_name.py' to start the visualization" << "\n";
        f_out << "#==========================================================================" << "\n" << "\n";
        
        bool lbp = false;
        bool pbp = false;
        bool cbp = false;
        bool wbp = false;
        bool mbp = false;
        bool lgp = false;
        bool pgp = false;
        bool cgp = false;
        bool wgp = false;
        bool mgp = false;
        
        //!Protein:
        float currpot,radius;
        int state = 1;
        int potindi[2] = {0,0};
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
            bool first = true;
            potindi[1] = 0;
            for (atoms_vec apro=grid_atoms.begin(); apro!=grid_atoms.end(); ++apro) {
                currpot = pro_vis.value(potindi);
                if (currpot > 0.) {
                    radius = currpot / vis_scaling;
                    if (radius < min_sphere_rad) radius = min_sphere_rad;
                    else if (radius > max_sphere_rad) radius = max_sphere_rad;
                    if (first) {f_out << "bad_pro_atom = [7.0,"; first = false;}
                    else f_out << ",7.0,";
                    if (flex_res) {
                        if (flex_map[state-1].find((*apro)->intern_id) == flex_map[state-1].end()) {
                            f_out << (*apro)->coord[0] << "," << (*apro)->coord[1] << "," 
                                  << (*apro)->coord[2] << "," << radius;
                        } else f_out << flex_map[state-1][(*apro)->intern_id][0] << ","
                                     << flex_map[state-1][(*apro)->intern_id][1] << ","
                                     << flex_map[state-1][(*apro)->intern_id][2] << "," << radius;
                    } else {
                        f_out << (*apro)->coord[0] << "," << (*apro)->coord[1] << "," 
                                  << (*apro)->coord[2] << "," << radius;
                    }
                }
                potindi[1] += 1;
            }
            if (!first) {
                f_out << "]" << "\n";
                f_out << "cmd.load_cgo(bad_pro_atom, 'pro_bad_potentials', " << state << ")" << "\n";
                pbp = true;
            } else {
                f_out << "bad_pro_atom = []\n"
                      << "cmd.load_cgo(bad_pro_atom, 'pro_bad_potentials', " << state << ")" << "\n";
                pbp = true;
            }
            ++state;
            potindi[0] += 1;
        }
        
        state = 1;
        potindi[0] = 0;
        potindi[1] = 0;
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
            bool first = true;
            potindi[1] = 0;
            for (atoms_vec apro=grid_atoms.begin(); apro!=grid_atoms.end(); ++apro) {
                currpot = pro_vis.value(potindi);
                if (currpot < 0.) {
                    radius = -currpot / vis_scaling;
                    if (radius < min_sphere_rad) radius = min_sphere_rad;
                    else if (radius > max_sphere_rad) radius = max_sphere_rad;
                    if (first) {f_out << "good_pro_atom = [7.0,"; first = false;}
                    else f_out << ",7.0,";
                    if (flex_res) {
                        if (flex_map[state-1].find((*apro)->intern_id) == flex_map[state-1].end()) {
                            f_out << (*apro)->coord[0] << "," << (*apro)->coord[1] << "," 
                                  << (*apro)->coord[2] << "," << radius;
                        } else f_out << flex_map[state-1][(*apro)->intern_id][0] << ","
                                     << flex_map[state-1][(*apro)->intern_id][1] << ","
                                     << flex_map[state-1][(*apro)->intern_id][2] << "," << radius;
                    } else {
                        f_out << (*apro)->coord[0] << "," << (*apro)->coord[1] << "," 
                                  << (*apro)->coord[2] << "," << radius;
                    }
                }
                potindi[1] += 1;
            }
            if (!first) {
                f_out << "]" << "\n";
                f_out << "cmd.load_cgo(good_pro_atom, 'pro_good_potentials', " << state << ")" << "\n";
                pgp = true;
            } else {
                f_out << "good_pro_atom = []\n";
                f_out << "cmd.load_cgo(good_pro_atom, 'pro_good_potentials', " << state << ")" << "\n";
                pgp = true;
            }
            ++state;
            potindi[0] += 1;
        }
        
        //!Ligand:
        state = 1;
        potindi[0] = 0;
        potindi[1] = 0;
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
            bool first = true;
            potindi[1] = 0;
            for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                currpot = lig_vis.value(potindi);
                if (currpot < 0.) {
                    radius = -currpot / vis_scaling;
                    if (radius < min_sphere_rad) radius = min_sphere_rad;
                    else if (radius > max_sphere_rad) radius = max_sphere_rad;
                    if (first) {f_out << "good_lig_atom = [7.0,"; first = false;}
                    else f_out << ",7.0,";
                    f_out << (*alig)->coord[0] << "," << (*alig)->coord[1] << "," << (*alig)->coord[2] << ","
                          << radius;
                }
                potindi[1] += 1;
            }
            if (!first) {
                f_out << "]" << "\n";
                f_out << "cmd.load_cgo(good_lig_atom, 'lig_good_potentials', " << state << ")" << "\n";
                lgp = true;
            } else {
                f_out << "good_lig_atom = []\n";
                f_out << "cmd.load_cgo(good_lig_atom, 'lig_good_potentials', " << state << ")" << "\n";
                lgp = true;
            }
            ++state;
            potindi[0] += 1;
        }
        
        state = 1;
        potindi[0] = 0;
        potindi[1] = 0;
        for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
            bool first = true;
            potindi[1] = 0;
            for (atoms_vec alig=(*lig)->atoms.begin(); alig!=(*lig)->atoms.end(); ++alig) {
                currpot = lig_vis.value(potindi);
                if (currpot > 0.) {
                    radius = currpot / vis_scaling;
                    if (radius < min_sphere_rad) radius = min_sphere_rad;
                    else if (radius > max_sphere_rad) radius = max_sphere_rad;
                    if (first) {f_out << "bad_lig_atom = [7.0,"; first = false;}
                    else f_out << ",7.0,";
                    f_out << (*alig)->coord[0] << "," << (*alig)->coord[1] << "," << (*alig)->coord[2] << ","
                          << radius;
                }
                potindi[1] += 1;
            }
            if (!first) {
                f_out << "]" << "\n";
                f_out << "cmd.load_cgo(bad_lig_atom, 'lig_bad_potentials', " << state << ")" << "\n";
                lbp = true;
            } else {
                f_out << "bad_lig_atom = []\n";
                f_out << "cmd.load_cgo(bad_lig_atom, 'lig_bad_potentials', " << state << ")" << "\n";
                lbp = true;
            }
            ++state;
            potindi[0] += 1;
        }
        
        //!Kofaktor:
        if (cof_file != "X") {
            state = 1;
            potindi[0] = 0;
            potindi[1] = 0;
            if (interaction_mode == 0 || interaction_mode == 5 ||
                interaction_mode == 6 || interaction_mode == 7) {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec acof=cof_structure->ligands[potindi[0]]->atoms.begin(); 
                                   acof!=cof_structure->ligands[potindi[0]]->atoms.end(); ++acof) {
                        currpot = cof_vis.value(potindi);
                        if (currpot < 0.) {
                            radius = -currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "good_cof_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*acof)->coord[0] << "," << (*acof)->coord[1] << "," << (*acof)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(good_cof_atom, 'cof_good_potentials', " << state << ")" << "\n";
                        cgp = true;
                    } else {
                        f_out << "good_cof_atom = []\n";
                        f_out << "cmd.load_cgo(good_cof_atom, 'cof_good_potentials', " << state << ")" << "\n";
                        cgp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            } else {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec acof=cof_structure->ligands[0]->atoms.begin(); 
                                   acof!=cof_structure->ligands[0]->atoms.end(); ++acof) {
                        currpot = cof_vis.value(potindi);
                        if (currpot < 0.) {
                            radius = -currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "good_cof_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*acof)->coord[0] << "," << (*acof)->coord[1] << "," << (*acof)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(good_cof_atom, 'cof_good_potentials', " << state << ")" << "\n";
                        cgp = true;
                    } else {
                        f_out << "good_cof_atom = []\n";
                        f_out << "cmd.load_cgo(good_cof_atom, 'cof_good_potentials', " << state << ")" << "\n";
                        cgp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            }
            state = 1;
            potindi[0] = 0;
            potindi[1] = 0;
            if (interaction_mode == 0 || interaction_mode == 5 ||
                interaction_mode == 6 || interaction_mode == 7) {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec acof=cof_structure->ligands[potindi[0]]->atoms.begin(); 
                                   acof!=cof_structure->ligands[potindi[0]]->atoms.end(); ++acof) {
                        currpot = cof_vis.value(potindi);
                        if (currpot > 0.) {
                            radius = currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "bad_cof_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*acof)->coord[0] << "," << (*acof)->coord[1] << "," << (*acof)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(bad_cof_atom, 'cof_bad_potentials', " << state << ")" << "\n";
                        cbp = true;
                    } else {
                        f_out << "bad_cof_atom = []\n";
                        f_out << "cmd.load_cgo(bad_cof_atom, 'cof_bad_potentials', " << state << ")" << "\n";
                        cbp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            } else {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec acof=cof_structure->ligands[0]->atoms.begin(); 
                                   acof!=cof_structure->ligands[0]->atoms.end(); ++acof) {
                        currpot = cof_vis.value(potindi);
                        if (currpot > 0.) {
                            radius = currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "bad_cof_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*acof)->coord[0] << "," << (*acof)->coord[1] << "," << (*acof)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(bad_cof_atom, 'cof_bad_potentials', " << state << ")" << "\n";
                        cbp = true;
                    } else {
                        f_out << "bad_cof_atom = []\n";
                        f_out << "cmd.load_cgo(bad_cof_atom, 'cof_bad_potentials', " << state << ")" << "\n";
                        cbp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            }
        }
        
        //!Metal:
        if (met_file != "X") {
            state = 1;
            potindi[0] = 0;
            potindi[1] = 0;
            if (interaction_mode == 0 || interaction_mode == 2 ||
                interaction_mode == 4 || interaction_mode == 6) {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec amet=met_structure->ligands[potindi[0]]->atoms.begin(); 
                                   amet!=met_structure->ligands[potindi[0]]->atoms.end(); ++amet) {
                        currpot = met_vis.value(potindi);
                        if (currpot < 0.) {
                            radius = -currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "good_met_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*amet)->coord[0] << "," << (*amet)->coord[1] << "," << (*amet)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(good_met_atom, 'met_good_potentials', " << state << ")" << "\n";
                        mgp = true;
                    } else {
                        f_out << "good_met_atom = []\n";
                        f_out << "cmd.load_cgo(good_met_atom, 'met_good_potentials', " << state << ")" << "\n";
                        mgp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            } else {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec amet=met_structure->ligands[0]->atoms.begin(); 
                                   amet!=met_structure->ligands[0]->atoms.end(); ++amet) {
                        currpot = met_vis.value(potindi);
                        if (currpot < 0.) {
                            radius = -currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "good_met_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*amet)->coord[0] << "," << (*amet)->coord[1] << "," << (*amet)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(good_met_atom, 'met_good_potentials', " << state << ")" << "\n";
                        mgp = true;
                    } else {
                        f_out << "good_met_atom = []\n";
                        f_out << "cmd.load_cgo(good_met_atom, 'met_good_potentials', " << state << ")" << "\n";
                        mgp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            }
            state = 1;
            potindi[0] = 0;
            potindi[1] = 0;
            if (interaction_mode == 0 || interaction_mode == 2 ||
                interaction_mode == 4 || interaction_mode == 6) {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec amet=met_structure->ligands[potindi[0]]->atoms.begin(); 
                                   amet!=met_structure->ligands[potindi[0]]->atoms.end(); ++amet) {
                        currpot = met_vis.value(potindi);
                        if (currpot > 0.) {
                            radius = currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "bad_met_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*amet)->coord[0] << "," << (*amet)->coord[1] << "," << (*amet)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(bad_met_atom, 'met_bad_potentials', " << state << ")" << "\n";
                        mbp = true;
                    } else {
                        f_out << "bad_met_atom = []\n";
                        f_out << "cmd.load_cgo(bad_met_atom, 'met_bad_potentials', " << state << ")" << "\n";
                        mbp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            } else {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec amet=met_structure->ligands[0]->atoms.begin(); 
                                   amet!=met_structure->ligands[0]->atoms.end(); ++amet) {
                        currpot = met_vis.value(potindi);
                        if (currpot > 0.) {
                            radius = currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "bad_met_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*amet)->coord[0] << "," << (*amet)->coord[1] << "," << (*amet)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(bad_met_atom, 'met_bad_potentials', " << state << ")" << "\n";
                        mbp = true;
                    } else {
                        f_out << "bad_met_atom = []\n";
                        f_out << "cmd.load_cgo(bad_met_atom, 'met_bad_potentials', " << state << ")" << "\n";
                        mbp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            }
        }
        
        //!Wasser:
        if (wat_file != "X" || gold_water) {
            state = 1;
            potindi[0] = 0;
            potindi[1] = 0;
            if (interaction_mode == 0 || interaction_mode == 3 ||
                interaction_mode == 4 || interaction_mode == 7 || gold_water) {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec awat=wat_structure->ligands[potindi[0]]->atoms.begin(); 
                                   awat!=wat_structure->ligands[potindi[0]]->atoms.end(); ++awat) {
                        currpot = wat_vis.value(potindi);
                        if (currpot < 0.) {
                            radius = -currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "good_wat_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*awat)->coord[0] << "," << (*awat)->coord[1] << "," << (*awat)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(good_wat_atom, 'wat_good_potentials', " << state << ")" << "\n";
                        wgp = true;
                    } else {
                        f_out << "good_wat_atom = []\n";
                        f_out << "cmd.load_cgo(good_wat_atom, 'wat_good_potentials', " << state << ")" << "\n";
                        wgp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            } else {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec awat=wat_structure->ligands[0]->atoms.begin(); 
                                   awat!=wat_structure->ligands[0]->atoms.end(); ++awat) {
                        currpot = wat_vis.value(potindi);
                        if (currpot < 0.) {
                            radius = -currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "good_wat_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*awat)->coord[0] << "," << (*awat)->coord[1] << "," << (*awat)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(good_wat_atom, 'wat_good_potentials', " << state << ")" << "\n";
                        wgp = true;
                    } else {
                        f_out << "good_wat_atom = []\n";
                        f_out << "cmd.load_cgo(good_wat_atom, 'wat_good_potentials', " << state << ")" << "\n";
                        wgp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            }
            state = 1;
            potindi[0] = 0;
            potindi[1] = 0;
            if (interaction_mode == 0 || interaction_mode == 3 ||
                interaction_mode == 4 || interaction_mode == 7 || gold_water) {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec awat=wat_structure->ligands[potindi[0]]->atoms.begin(); 
                                   awat!=wat_structure->ligands[potindi[0]]->atoms.end(); ++awat) {
                        currpot = wat_vis.value(potindi);
                        if (currpot > 0.) {
                            radius = currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "bad_wat_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*awat)->coord[0] << "," << (*awat)->coord[1] << "," << (*awat)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(bad_wat_atom, 'wat_bad_potentials', " << state << ")" << "\n";
                        wbp = true;
                    } else {
                        f_out << "bad_wat_atom = []\n";
                        f_out << "cmd.load_cgo(bad_wat_atom, 'wat_bad_potentials', " << state << ")" << "\n";
                        wbp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            } else {
                for (ligands_vec lig=lig_structure->ligands.begin(); lig!=lig_structure->ligands.end();++lig) {
                    bool first = true;
                    potindi[1] = 0;
                    for (atoms_vec awat=wat_structure->ligands[0]->atoms.begin(); 
                                   awat!=wat_structure->ligands[0]->atoms.end(); ++awat) {
                        currpot = wat_vis.value(potindi);
                        if (currpot > 0.) {
                            radius = currpot / vis_scaling;
                            if (radius < min_sphere_rad) radius = min_sphere_rad;
                            else if (radius > max_sphere_rad) radius = max_sphere_rad;
                            if (first) {f_out << "bad_wat_atom = [7.0,"; first = false;}
                            else f_out << ",7.0,";
                            f_out << (*awat)->coord[0] << "," << (*awat)->coord[1] << "," << (*awat)->coord[2] << ","
                            << radius;
                        }
                        potindi[1] += 1;
                    }
                    if (!first) {
                        f_out << "]" << "\n";
                        f_out << "cmd.load_cgo(bad_wat_atom, 'wat_bad_potentials', " << state << ")" << "\n";
                        wbp = true;
                    } else {
                        f_out << "bad_wat_atom = []\n";
                        f_out << "cmd.load_cgo(bad_wat_atom, 'wat_bad_potentials', " << state << ")" << "\n";
                        wbp = true;
                    }
                    ++state;
                    potindi[0] += 1;
                }
            }
        }
        
        //!Jetzt die good- und bad-Distances:
        for (map<int,vector<TRIPLE> >::iterator mit=good_con.begin(); mit!=good_con.end(); ++mit) {
            f_out << "good_inter =  [";
            for (vector<TRIPLE>::iterator mjt=mit->second.begin(); mjt!=mit->second.end(); ++mjt) {
                if (mjt == mit->second.begin()) f_out << "9.0,";
                else f_out << ",9.0,";
                f_out << mjt->a[0] << "," << mjt->a[1] << "," << mjt->a[2] << ","
                      << mjt->b[0] << "," << mjt->b[1] << "," << mjt->b[2] << ","
                      << crad << ",0.0,0.0,1.0,0.0,0.0,1.0";
            }
            f_out << "]" << "\n";
            f_out << "cmd.load_cgo(good_inter, 'good_distances', " << mit->first + 1 << ")" << "\n";
        }
        for (map<int,vector<TRIPLE> >::iterator mit=bad_con.begin(); mit!=bad_con.end(); ++mit) {
            f_out << "bad_inter =  [";
            for (vector<TRIPLE>::iterator mjt=mit->second.begin(); mjt!=mit->second.end(); ++mjt) {
                if (mjt == mit->second.begin()) f_out << "9.0,";
                else f_out << ",9.0,";
                f_out << mjt->a[0] << "," << mjt->a[1] << "," << mjt->a[2] << ","
                      << mjt->b[0] << "," << mjt->b[1] << "," << mjt->b[2] << ","
                      << crad << ",1.0,0.0,0.0,1.0,0.0,0.0";
            }
            f_out << "]" << "\n";
            f_out << "cmd.load_cgo(bad_inter, 'bad_distances', " << mit->first + 1 << ")" << "\n";
        }
        
        //Farben:
        if (lbp) f_out << "cmd.color('red', 'lig_bad_potentials')" << "\n";
        if (pbp) f_out << "cmd.color('red', 'pro_bad_potentials')" << "\n";
        if (cof_file != "X" && cbp) f_out << "cmd.color('red', 'cof_bad_potentials')" << "\n";
        if (met_file != "X" && mbp) f_out << "cmd.color('red', 'met_bad_potentials')" << "\n";
        if ((wat_file != "X" || gold_water) && wbp) f_out << "cmd.color('red', 'wat_bad_potentials')" << "\n";
        if (lgp) f_out << "cmd.color('blue', 'lig_good_potentials')" << "\n";
        if (pgp) f_out << "cmd.color('blue', 'pro_good_potentials')" << "\n";
        if (cof_file != "X" && cgp) f_out << "cmd.color('blue', 'cof_good_potentials')" << "\n";
        if (met_file != "X" && mgp) f_out << "cmd.color('blue', 'met_good_potentials')" << "\n";
        if ((wat_file != "X" || gold_water) && wgp) f_out << "cmd.color('blue', 'wat_good_potentials')" << "\n";
        
        //Protein und Liganden laden:
        f_out << "#--------------------- loading protein and ligand -------------------------" << "\n";
        if (flex_res) f_out << "cmd.load('" << fname << "','protein')" << "\n";
        else f_out << "cmd.load('" << pro_file << "','protein')" << "\n";
        if (lig_type == "dlg") {
            f_out << "cmd.load('" << newmolname << "','ligands')" << "\n";
        } else f_out << "cmd.load('" << lig_file << "','ligands')" << "\n";
        if (cof_file != "X") {
            f_out << "cmd.load('" << cof_file << "','cofactor')" << "\n";
            f_out << "cmd.show('sticks','cofactor')" << "\n";
        }
        
        if (met_file != "X") {
            f_out << "cmd.load('" << met_file << "','metals')" << "\n";
            f_out << "cmd.show('lines','metals')" << "\n";
        }
        if (wat_file != "X"/* || gold_water*/) {
            f_out << "cmd.load('" << wat_file << "','waters')" << "\n";
            f_out << "cmd.show('lines','waters')" << "\n";
        } else if (wname != "X") {
            f_out << "cmd.load('" << wname << "','waters')" << "\n";
            f_out << "cmd.show('lines','waters')" << "\n";
        }
        
        f_out << "cmd.set('stick_radius','0.06')" << "\n";
        f_out << "cmd.show('sticks','ligands')" << "\n";
        f_out << "cmd.select('pocket','ligands expand " << flex_cav_rad << "')" << "\n";
        f_out << "cmd.show('sticks','pocket')" << "\n";
        if (!flex_res) f_out << "cmd.select('far','protein AND NOT pocket')" << "\n";
        if (!flex_res) f_out << "cmd.hide('everything', 'far')" << "\n";
        f_out << "cmd.select('none')" << "\n";
        f_out << "cmd.zoom('pocket')" << "\n";
        f_out << "cmd.hide('everything','elem LP')" << "\n";
        f_out << "cmd.hide('everything','elem h')" << endl;
        
        if (score_torsions) {
            for (map<int,vector<int> >::iterator tors=bad_tors.begin(); tors!=bad_tors.end(); ++tors) {
                int n_tors = tors->second.size() / 4;
                for (int nt=0; nt<n_tors; ++nt) {
                //    f_out << "cmd.set('state'," << tors->first+1 << ")" << "\n";
                    f_out << "cmd.select('pk1','id " << tors->second[4*nt] << " AND ligands')" << "\n";
                    f_out << "cmd.select('pk2','id " << tors->second[4*nt+1] << " AND ligands')" << "\n";
                    f_out << "cmd.select('pk3','id " << tors->second[4*nt+2] << " AND ligands')" << "\n";
                    f_out << "cmd.select('pk4','id " << tors->second[4*nt+3] << " AND ligands')" << "\n";
                    f_out << "cmd.dihedral('bad_torsions',state='" << tors->first+1 << "')" << "\n";
                }
            }
        }
        
        f_out.close();
    } else if (visualization || (n_visi < 55558)) {
        if (visualization) {
            if (!silent_mode) cout << " --> calculating visualization elements and writing files for top " << n_visi << endl;
        }
        string topname = lig_file;
        multimap<float,int> ssorted; //score : number   der gleiche Score kann mehrfach vorkommen => Multimap
        for (map<int,float>::iterator it=score_map.begin(); it!=score_map.end(); ++it) {
            ssorted.insert(make_pair(it->second,it->first));
        }
        stl_ptr<STRUCTURE> tmp_structure;
        stl_ptr<PARSER> tmp_parser;
        tmp_structure = new STRUCTURE(0);
        tmp_parser = new PARSER(&(*tmp_structure),0,false,true,true,false);
        int zz = 0;
        for (multimap<float,int>::iterator it=ssorted.begin(); it!=ssorted.end(); ++it) {
            stl_ptr<LIGAND> plig;
            plig = new LIGAND(*(lig_structure->ligands[it->second]));
            plig->atom_types_already_set = false;
            plig->get_atom_typing(1,false,"X",true,0,true,false,"X");
            tmp_structure->ligands.push_back(plig);
            ++zz;
            if (zz == n_visi) break;
        }
        
        string newext = "_top";
        string_fu::add2s(newext,n_visi);
        string_fu::add2s(newext,".mol2");
        string_fu::replace_ext(topname,".mol2",newext);
        if (!visualization) {
            if (!silent_mode) cout << " --> writing top " << n_visi << " ligands file: " << topname << endl;
        }
        tmp_parser->write_mol2(topname.c_str());
        tmp_structure.kill();
        tmp_parser.kill();
    } else {
        string wname = "X";
        //!neu 29.07.2010:  Tobi will die GOLD-water auch als mol2 rausgeschrieben haben: ==========================
        if (gold_water) {
            vector<stl_ptr<LIGAND> > ltw;
            for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                    (*at)->sybyl_type = "O.3";
                }
                ltw.push_back(*lig);
            }
            wname = lig_file;
            string_fu::replace_ext(wname,".mol2","_gold_water.mol2");
            if (!silent_mode) cout << " --> writing " << wname << " for visualization" << endl;
            lig_parser->write_mol2(ltw,wname.c_str());
        }
        //!=========================================================================================================
    }
}
//!============================================================================================================================


//!============================================================================================================================
void show_header() {
    cout << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n";
    cout << " | 'DSX'           Knowledge-based scoring function for the assessment       |" << "\n";
    cout << " |                 of receptor-ligand interactions                           |" << "\n";
    cout << " |  author     :   Gerd Neudert                                              |" << "\n";
    cout << " |  supervisor :   Prof. Dr. G. Klebe                       " << c_message<cGREY>("___    _ _") << "       |" << "\n";
    cout << " |  mailto     :   gerd.neudert@gmail.com                   " << c_message<cGREY>("))_    )\\`)") << "      |" << "\n";
    cout << " |  version    :   " << version_string << "   (" << version_date << ")                     " << c_message<cGREY>("((_( o ((\\( o") << "     |" << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n" << "\n";
    cout << copyright_string;
    cout << " => Call:" << "\n";
    cout << "  dsx [-cghjropsv] -P pro_file -L lig_file [-I int] [-S int] [-Tx int]\n"
         << "      [-Vx float] [-C cof_file]  [-W wat_file] [-M met_file]\n"
         << "      [-R ref_file]  -D pot_dir\n";
    cout << "========================================================================" << "\n" << "\n";
}

void show_help() {
    show_header();
    cout << "Type 'dsx -h' to get more help" << "\n" << endl;
}

void show_big_help() {
    show_header();
    cout << "  pro_file    :  A pdb or mol2 file of your protein.\n";
    cout << "                  In pdb format metals in this file will be treated as part\n"
         << "                  of the protein. => Be sure to delete metals in the pdb file\n"
         << "                  if you want to supply some metals seperately (-M met_file)!\n";
    cout << "                  All other HETATMs will be ignored!\n";
    cout << "                  In mol2 format everything will be taken as part of the\n"
         << "                  protein. => Be sure to delete molecules you want to supply\n"
         << "                  seperately (-C, -W, -M) from the protein-mol2-file!\n";
    cout << "  lig_file    :  A mol2- or autodock dlg-file containing all molecules that\n"
         << "                  should be scored.\n";
    cout << "  cof_file    :  A mol2 file supplying a cofactor that should be used.\n";
    cout << "  wat_file    :  A multimol2 file supplying all water molecules that should\n"
         << "                  be considered.\n";
    cout << "  met_file    :  A multimol2 file supplying all metal atoms that should\n"
         << "                  be considered.\n";
    cout << "  ref_file    :  A mol2 file used as reference for rmsd calculation.\n";
    cout << "  pot_dir     :  The directory with the potentials to use.\n";
    cout << "                  Alternativly you can set the DSX_POTENTIALS environment\n"
         << "                  variable to the pot_dir.\n";
    cout << "  -c          :  Turns ON the check for covalently bound ligands. Ligands\n"
         << "                  sharing atom positions with the protein will be considered\n"
         << "                  as covalently bound ligands (via the shared atoms). Please\n"
         << "                  note that the scoring regarding covalent bond ligands\n"
         << "                  is much slower.\n";
    cout << "  -f          :  Scoring of docking solutions with flexible protein parts.\n"
         << "                  If you performed a flexible docking with AUTODOCK you have\n"
         << "                  to supply the corresponding dlg-file instead of a mol2.\n";
    cout << "                  For GOLD solutions the necessary data should be inside\n"
         << "                  the mol2. This mode is currently only available with\n"
         << "                  pdb-proteins.\n";
    cout << "  -g          :  Use this flag if you want to score docking solutions from\n"
         << "                  a GOLD run considering water molecules. For each single\n"
         << "                  solution the water molecules marked as 'ON' will be scored.\n";
    cout << "  -h          :  Prints this help.\n";
    cout << "  -j          :  Minimize Ligands to a close local minimum.\n";
    cout << "                  If you use this option you should also use the '-t' flag\n"
         << "                  for better results.\n";
    cout << "  -r          :  Turns on rmsd calculation. If you don't specify a reference\n"
         << "                  structure the first molecule from your lig_file will be used\n"
         << "                  as reference. Values of '-1' indicate a failure in rmsd\n"
         << "                  calculation, e.g. if two different molecules are compared.\n";
    cout << "  -o          :  Use this flag as an alternative rmsd calculation." << "\n";
    cout << "                  Normally in '-r' mode only equal molecules can be compared.\n"
         << "                  Due to disconnected bonds or other reasons it is possible\n"
         << "                  that some geometries become unequal compared to the\n"
         << "                  reference. In such cases '-r' would return '-1', but in\n"
         << "                  '-o' mode the maximum common substructure is matched. (only\n"
         << "                  in those cases where the normal function retruns -1).\n";
    cout << "  -p          :  Writes out aminoacid-profiles. (experimental)\n";
    cout << "  -s          :  Turns on the silent mode (output only errors).\n";
    cout << "  -v          :  Writes out a pymol visualization file.\n";
    cout << "                  A python file, which can be run from pymol and will show\n"
         << "                  you blue and red spheres for 'good' and 'bad' potentials\n"
         << "                  on each atom. Additionally it will show you blue and red\n"
         << "                  lines indicating single atom-atom interactions with a\n"
         << "                  potential lower or higher than a given threshold.\n";
    cout << "  -I int      :  Here you can specify the mode that affects how cofactors,\n"
         << "                  waters and metals will be handeled.\n";
    cout << "                  The default mode is '-I 1', which means, that all molecules\n"
         << "                  are treated as part of the protein. If a structure should\n"
         << "                  not be treated as part of the protein you have supply a\n"
         << "                  seperate file with seperate MOLECULE entries corresponding\n"
         << "                  to each MOLECULE entry in the ligand_file (It is assumed\n"
         << "                  that the structure, e.g. a cofactor, was kept flexible in\n"
         << "                  docking, so that there should be a different geometry\n"
         << "                  corresponding to each solution. Otherwise it won't make\n"
         << "                  sense not to treat it as part of the protein.).\n";
    cout << "                  The following modes are possible: " << "\n";
    cout << "                    0: cofactors, waters and metals interact with protein,\n"
         << "                       ligand and each other\n";
    cout << "                    1: cofactors, waters and metals are treated as part of\n"
         << "                       the protein\n";
    cout << "                    2: cofactors and metals are treated as part of the protein\n"
         << "                       (waters as in mode 0)\n";
    cout << "                    3: cofactors and waters are treated as part of the protein\n";
    cout << "                    4: cofactors are treated as part of the protein\n";
    cout << "                    5: metals and waters are treated as part of the protein\n";
    cout << "                    6: metals are treated as part of the protein\n";
    cout << "                    7: waters are treated as part of the protein\n";
    cout << "                  Please note: Only those structures can be treated\n"
         << "                  individually, which are supplied in seperate files.\n";
    cout << "  -S int      :  Here you can specify the mode that affects how the results\n"
         << "                  will be sorted. The default mode is '-S 1', which sorts the\n"
         << "                  ligands in the same order as they are found in the lig_file.\n"
         << "                  The following modes are possible:\n";
    cout << "                    0: Same order as in the ligand file\n";
    cout << "                    1: Ordered by increasing total score\n";
    cout << "                    2: Ordered by increasing per-atom-score\n";
    cout << "                    3: Ordered by increasing per-contact-score\n";
    cout << "                    4: Ordered by increasing rmsd\n";
    cout << "                    5: Ordered by increasing torsion score\n";
    cout << "                    6: Ordered by increasing per-torsion-score\n";
    cout << "  -Tx float   :  Specify the weightings for the different scoring terms.\n";
    cout << "                  The integer x specifies the exact term. Possible number\n";
    cout << "                  and corresponding default weights are listed below.\n";
    cout << "                    T0: distance dependent pair potentials (default: 1.0)\n";
    cout << "                    T1: intramolecular clashes             (default: 0.0)\n";
    cout << "                    T2: torsion potentials                 (default: 0.0)\n";
    cout << "                    T3: sas potentials                     (default: 0.0)\n";
    cout << "  -V1 float   :  Specify the scaling factor for the spheres in pymol\n"
         << "                  visualization (default: 16.). Increasing that number will\n"
         << "                  decrease the sphere radius.\n";
    cout << "  -V2 float   :  Specify the lower threshold for atom-atom interactions that\n"
         << "                  should be visualized (default: -1.).\n";
    cout << "  -V3 float   :  Specify the upper threshold for atom-atom interactions that\n"
         << "                  should be visualized (default: 1.).\n";
    cout << "  -V4 int     :  If you don't want a visualization for all molecules, here you\n"
         << "                  can specify the number of top-ranked molecules that should\n"
         << "                  be visualized.\n";
    cout << "  -V5 float   :  Specify the threshold for bad torsion angles that should be\n"
         << "                  visualized (default: 1.).\n" << endl;
}

int main(int argc,char *argv[]) {
    
    vector<string> s;
    s.resize(argc);
    for (int i=1; i<argc; ++i) { //Kommandozeilenparameter in strings umwandeln
        s[i].assign(argv[i]);
    }
    for (int i=1; i<argc; ++i) { //Leerzeichen entfernen
        ostringstream os;
        istringstream is;
        is.clear(); is.str(s[i]);
        string tmp; is >> tmp;
        os << tmp; s[i] = os.str();
    }
    
    //!---------------------------------------------------------------
    //!Parameter:
    string pro_file = "X"; //Name des Protein Files (pdb oder mol2)
    string lig_file = "X"; //Name des Ligand Files (mol2 oder dlg)
    string cof_file = "X"; //Name des optionalen Kofaktor Files (mol2)
    string wat_file = "X"; //Name des optionalen Wasser Files (mol2)
    string met_file = "X"; //Name des optionalen Metall Files (mol2)
    string ref_file = "X"; //Name des optionalen Referenzstruktur Files (mol2)
    string pot_dir = "X"; //Name des Verzeichnis, in dem sie die zu benutzenden Potentiale befinden
    string pro_type = "pdb"; //Dateityp des Protein Files
    string lig_type = "mol"; //Dateityp des Ligandfile
    
    bool coval_check = false; //Auf kovalent gebundene Liganden pruefen?
    bool pairs; // Paarpotentiale nutzen
    bool torsions; // Torsionswinkel scoren?
    bool intra; // Intramolekulare Clashes scoren?
    bool sas; // SAS-Term berechnen
    bool gold_water = false; //Wassermolekuele aus einem GOLD-Docking beruecksichtigen?
    bool jiggle = false; //locale Minimierung der Liganden?
    bool jiggle_first_only = false; // Debug-Flag (undokumentiert): Nur die erste Pose Minimieren
    bool calc_rmsd = false; //RMSD-Werte berechnen?
    bool mcis = false; //RMSD mit Bron-Kerbosch berechnen, wenn normaler Mode -1 liefert
    bool silent_mode = false; //Nur Errors und Warnings ausgeben?
    bool visualization = false; //Python Files zur Visualisierung in Pymol rausschreiben?
    bool flex_res = false; //Verschiedene Proteinkonformationen beruecksichtigen?
    bool as_profile = false; //Sollen Aminosaeureprofile rausgeschrieben werden?
    
    int interaction_mode = 1; //Festlegen welche Atom-Atom Paarungen moeglich sind
    int sort_mode = 1; //In welcher Reihenfolge sollen die Ergebnisse ausgegeben werden
    int n_visi = 99999; //Anzahl der Liganden, fuer die eine Visualisierung rausgeschrieben werden soll
    
    float vis_scaling = 16.; //Skalierungsfaktor fuer die roten und blauen Kugeln
    float vis_lower = -1.; //Atom-Atom Interaktionen mit einem Potential unter dieser Schwelle visualisieren
    float vis_upper = 1.; //Atom-Atom Interaktionen mit einem Potential oberhalb dieser Schwelle visualisieren
    float vis_tors = 1.;

    //!---------------------------------------------------------------
    
    //!---------------------------------------------------------------
    //!Kommandozeile parsen:
    for (int i=1; i<argc; ++i) {
        if (s[i][0] == '-') {
            for (int j=1; j<int(s[i].size()); ++j) {
                if (s[i][j] == 'c') coval_check = true;
                else if (s[i][j] == 'f') flex_res = true;
                else if (s[i][j] == 'g') gold_water = true;
                else if (s[i][j] == 'j') jiggle = true;
                else if (s[i][j] == 'k') jiggle_first_only = true;
                else if (s[i][j] == 'r') calc_rmsd = true;
                else if (s[i][j] == 'o') {calc_rmsd = true; mcis = true;}
                else if (s[i][j] == 'p') as_profile = true;
                else if (s[i][j] == 's') silent_mode = true;
                else if (s[i][j] == 'v') visualization = true;
                else if (s[i][j] == 'h') {
                    show_big_help();
                    exit(0);
                }
                else if (s[i][j] == 'I') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -I") << endl;
                        show_help();
                        exit(1);
                    } else {
                        if (!string_fu::s2v(s[i],interaction_mode)) {
                            cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                            show_help();
                            exit(1);
                        }
                        if (interaction_mode < 0 || interaction_mode > 7) {
                            cerr << c_message<cWARNING>(s[i]) << " is not a valid mode! -> using default instead" << endl;
                            interaction_mode = 0;
                        }
                    }
                    break;
                }
                else if (s[i][j] == 'S') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -S") << endl;
                        show_help();
                        exit(1);
                    } else {
                        if (!string_fu::s2v(s[i],sort_mode)) {
                            cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                            show_help();
                            exit(1);
                        }
                        if (sort_mode < 0 || sort_mode > 6) {
                            cerr << c_message<cWARNING>(s[i]) << " is not a valid mode! -> using default instead" << endl;
                            sort_mode = 0;
                        }
                    }
                    break;
                }
                else if (s[i][j] == 'T') {
                    if ((j+1) == int(s[i].size())) {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (use T0, T1, T2 or T3)!" << endl;
                        show_help();
                        exit(1);
                    }
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -Tx") << endl;
                        show_help();
                        exit(1);
                    } else {
                        if (s[i-1][j+1] == '0') {
                            if (!string_fu::s2v(s[i],pair_weight)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '1') {
                            if (!string_fu::s2v(s[i],intra_weight)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '2') {
                            if (!string_fu::s2v(s[i],torsion_weight)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '3') {
                            if (!string_fu::s2v(s[i],sas_weight)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        }
                    }
                    break;
                }
                else if (s[i][j] == 'V') {
                    if ((j+1) == int(s[i].size())) {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (use V1, V2, ... or V5)!" << endl;
                        show_help();
                        exit(1);
                    }
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -Vx") << endl;
                        show_help();
                        exit(1);
                    } else {
                        if (s[i-1][j+1] == '1') {
                            if (!string_fu::s2v(s[i],vis_scaling)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '2') {
                            if (!string_fu::s2v(s[i],vis_lower)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '3') {
                            if (!string_fu::s2v(s[i],vis_upper)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '4') {
                            if (!string_fu::s2v(s[i],n_visi)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        } else if (s[i-1][j+1] == '5') {
                            if (!string_fu::s2v(s[i],vis_tors)) {
                                cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                                show_help();
                                exit(1);
                            }
                        }
                    }
                    break;
                }
                else if (s[i][j] == 'P') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -P") << endl;
                        show_help();
                        exit(1);
                    } else {
                        pro_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'L') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -L") << endl;
                        show_help();
                        exit(1);
                    } else {
                        lig_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'C') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -C") << endl;
                        show_help();
                        exit(1);
                    } else {
                        cof_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'W') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -W") << endl;
                        show_help();
                        exit(1);
                    } else {
                        wat_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'M') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -M") << endl;
                        show_help();
                        exit(1);
                    } else {
                        met_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'R') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -R") << endl;
                        show_help();
                        exit(1);
                    } else {
                        ref_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'D') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -D") << endl;
                        show_help();
                        exit(1);
                    } else {
                        pot_dir = s[i];
                    }
                    break;
                }
                else {
                    cerr << c_message<cERROR>("unknown option -") << s[i][j] << endl;
                    show_help();
                    exit(1);
                }
            }
        }
    }
    //!---------------------------------------------------------------
    
    if (argc < 5) {
        cerr << c_message<cERROR>("not enough arguments") << endl;
        show_help();
        exit(1);
    }
    
    //!---------------------------------------------------------------
    //!Gucken, ob mol2- oder pdb-Protein:
    if (pro_file == "X") {
        cerr << c_message<cERROR>("you have to specify a protein (pdb or mol2)") << endl;
        show_help();
        exit(1);
    }
    string pro_ext = string_fu::get_ext(pro_file);
    if (pro_ext == ".mol2" || pro_ext == ".MOL2") {
        pro_type = "mol";
    } else pro_type = "pdb";
    //!---------------------------------------------------------------
    
    //!---------------------------------------------------------------
    //!Gucken, ob mol2- oder dlg-Ligandenfile:
    if (lig_file == "X") {
        cerr << c_message<cERROR>("you have to specify a ligands file (multimol2)") << endl;
        show_help();
        exit(1);
    }
    string lig_ext = string_fu::get_ext(lig_file);
    if (lig_ext == ".mol2" || lig_ext == ".MOL2") lig_type = "mol";
    else lig_type = "dlg";
    //!---------------------------------------------------------------
    
    //!---------------------------------------------------------------
    //!Check, ob die Kommandozeile valide ist:
    if (pot_dir == "X") {
        char *envbuf;
        envbuf = getenv("DSX_POTENTIALS");
        if (envbuf != NULL) {
            pot_dir = envbuf;
        } else {
            cerr << c_message<cERROR>("you have to specify the directory with the set of potentials to use") << endl;
            cerr << "       (alternativly you can set the DSX_POTENTIALS environment variable to point on it)" << endl;
            show_help();
            exit(1);
        }
    }
    
    if (flex_res) {
        if (pro_type != "pdb") {
            cerr << c_message<cERROR>("flexible residues are currently only possible with protein in pdb format") << endl;
            show_help();
            exit(1);
        }
    }
    //!---------------------------------------------------------------

    if (pair_weight > weight_limit) pairs = true;
    else pairs = false;
    if (intra_weight > weight_limit) intra = true;
    else intra = false;
    if (torsion_weight > weight_limit) torsions = true;
    else torsions = false;
    if (sas_weight > weight_limit) sas = true;
    else sas = false;

    
    if (((cof_file != "X" || wat_file != "X" || met_file != "X") && pro_type == "mol") ||
        ((met_file != "X") && pro_type == "pdb")) {
        if (!silent_mode) {
            cout << "\n" << c_message<cWARNING>("you supplied additional input files. Please make sure that the supplied structures")
                 << "\n";
            cout << "         are not also included in your protein file!!!" << endl;
        }
    } else if (gold_water && pro_type == "mol") {
        if (!silent_mode) {
            cout << "\n" << c_message<cWARNING>("reading GOLD water is switched on and your protein is a mol2 file. ")
                 << "Please make sure to have no waters in the protein file!" << endl;
        }
    }
    
    if (jiggle && coval_check) {
        cerr << c_message<cERROR>("minimization of covalently bond ligands is currently not possible.") << "\n";
        cerr << "Please make sure to have no covalent ligands and restart without '-c'." << endl;
        exit(1);
    }
    
    if (!silent_mode) cout << endl;
    
    //!---------------------------------------------------------------
    //!Scoring-Objekt initialisieren:
    stl_ptr<SCORER> myscore = new SCORER(pro_file,lig_file,cof_file,wat_file,met_file,ref_file,pot_dir,pro_type,lig_type,
                                         coval_check,pairs,torsions,intra,sas,gold_water,jiggle,jiggle_first_only,calc_rmsd,mcis,silent_mode,visualization,
                                         flex_res,as_profile,interaction_mode,sort_mode,vis_scaling,vis_lower,vis_upper,vis_tors);
    //!---------------------------------------------------------------
    
    //!---------------------------------------------------------------
    //!Scoren:
    myscore->score();
    //!---------------------------------------------------------------
    
    //!---------------------------------------------------------------
    //!Ergebnisse schreiben:
    myscore->write_results(n_visi); //Visualisierung hier mit drinn
    if ((n_visi < 99999) && visualization) {
    
        if (jiggle) {
            string jname = lig_file;
            string_fu::replace_ext(jname,".mol2","_jiggled.mol2");
            lig_file = jname;
            jiggle = false;
        }
        
        myscore.kill();
        string topname = lig_file;
        string newext = "_top";
        string_fu::add2s(newext,n_visi);
        string_fu::add2s(newext,".mol2");
        string_fu::replace_ext(topname,".mol2",newext);
        myscore = new SCORER(pro_file,topname,cof_file,wat_file,met_file,ref_file,pot_dir,pro_type,lig_type,
                                 coval_check,pairs,torsions,intra,sas,gold_water,jiggle,jiggle_first_only,calc_rmsd,mcis,true,true,flex_res,as_profile,interaction_mode,
                                 sort_mode,vis_scaling,vis_lower,vis_upper,vis_tors);
        myscore->score();
        myscore->write_results(99999,false);
    }
    //!---------------------------------------------------------------
    
    //!---------------------------------------------------------------
    //!Aufraeumen:
    myscore.kill();
    //!---------------------------------------------------------------
    
    if (!silent_mode) cout << "done!" << endl;
    
    return 0;
}


