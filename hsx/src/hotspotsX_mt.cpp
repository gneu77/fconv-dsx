// Copyright (c)  Gerd Neudert 2007, 2008, 2009, 2011

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include"structure_GN.h"
#include"files_GN.h"
#include"linalg_GN.hpp"
#include"grid_GN.hpp"
#include"stl_ptr_GN.hpp"
#include<vector>
#include<tr1/unordered_map>
#include<algorithm>
#include<pthread.h>
#include<semaphore.h>
#include<stack>

#include"message_GN.hpp"
#include"string_fu_GN.hpp"
#include"gauss_GN.hpp"

#include"atom_GN.h"
#include"molecule_GN.h"
#include"protein_GN.h"
#include"structure_additional_GN.h"
#include"atom_properties_GN.h"

#if defined (_SPEEDTEST)
#include"timer_GN.hpp"
TIMER_GN<double> t_total;
#endif

using namespace std;
using namespace TXT;

const string copyright_string =
" Copyright (C) 2009, 2010, 2011 Gerd Neudert and Gerhard Klebe\n\
  Usage of HotspotsX is free without any limitations.\n\
  There is NO warranty, not even for MERCHANTABILITY or FITNESS\n\
  FOR A PARTICULAR PURPOSE.\n\n";


namespace MT {

    #define USE_N_THREADS 4

    pthread_attr_t attr;
    int smul1;
    int smul2;
    int ind1;
    int ind2;
    int ind3;
    tr1::unordered_map<string,vector<string> > rel_key_map;
    vector<stl_ptr<ATOM> > grid_atoms;
    vec3d<float> gmin,gmax;
    int gsize[3];
    int gs;
    vec3d<float>* acnt_grid;
    float max_square_dst;
    //GRID<2,float> sgrid;
    int s0;
    float* gsgrid;
    tr1::unordered_map<string,int> key_map;

    pthread_t wt[USE_N_THREADS];
    float spacing;
    float max_dst;
    bool normalize;
    bool good_bad;
    bool smooth;
    string acnt_mode;

    sem_t work_signal;
    stack<string> key_stack;
    sem_t master_wait;
    sem_t total_keys_worked;
    pthread_mutex_t wait_mutex = PTHREAD_MUTEX_INITIALIZER;
    int keys_done;
    int keys_todo;

    void* key_box(void* arg) {
        while (true) {
            sem_wait(&work_signal);

            string lkey = key_stack.top();
            key_stack.pop();

            bool last_one = false;
            if (key_stack.size() > 0) {
                sem_post(&work_signal);
            } else last_one = true;

            cout << lkey << "  "; cout.flush();
            int grs[3] = {gsize[0],gsize[1],gsize[2]};
            float *sg = new float[gs]; memset(sg,0,gs*sizeof(float));
            int g_ind2[3];
            float min_score = 0.;
            float max_score = 0.;
            int from_z,to_z,from_y,to_y,from_x,to_x;
            int steps = int((max_dst / spacing) + 1.01);
            
            for (vector<string>::iterator vkey=rel_key_map[lkey].begin(); vkey!=rel_key_map[lkey].end(); ++vkey) {
                //!Ueber alle Proteinatome:
                for (atoms_vec apro=grid_atoms.begin(); apro!=grid_atoms.end(); ++apro) {
                    if ((*apro)->sybyl_type == "X") continue;
                    string ckey = *vkey + "_" + (*apro)->sybyl_type;
                    if (key_map.find(ckey) == key_map.end()) continue;
                    
                    //!Den naechstgelegenen Gitterpunkt vom acnt_grid ermitteln:
                    g_ind2[0] = int((((*apro)->coord[0] - gmin[0]) / spacing) + 0.5);
                    g_ind2[1] = int((((*apro)->coord[1] - gmin[1]) / spacing) + 0.5);
                    g_ind2[2] = int((((*apro)->coord[2] - gmin[2]) / spacing) + 0.5);

                    from_z = g_ind2[2] - steps; to_z = g_ind2[2] + steps + 1;
                    if (from_z < 0) from_z = 0; if (to_z > grs[2]) to_z = grs[2];
                    from_y = g_ind2[1] - steps; to_y = g_ind2[1] + steps + 1;
                    if (from_y < 0) from_y = 0; if (to_y > grs[1]) to_y = grs[1];
                    from_x = g_ind2[0] - steps; to_x = g_ind2[0] + steps + 1;
                    if (from_x < 0) from_x = 0; if (to_x > grs[0]) to_x = grs[0];
                    //!Ueber die Box laufen:
                    for (int z=from_z; z<to_z; ++z) {
                        for (int y=from_y; y<to_y; ++y) {
                            for (int x=from_x; x<to_x; ++x) {
                                int gind = x; gind += y*smul1; gind += z*smul2;
                                
                                float sdst = get_square_distance(acnt_grid[gind],(*apro)->coord);
                                if (sdst > max_square_dst) continue;

                                sdst = sqrt(sdst);
                                int bin_key = int(sdst*100.);

                                int sind = key_map[ckey]; sind += bin_key*s0;
                                sg[gind] += gsgrid[sind];
                            }
                        }
                    }
                }
            }

            if (smooth) {
                float *fsg = new float[gs]; memset(fsg,0,gs*sizeof(float));

                float sigma = 0.8 * spacing;
                float gauss_orig = 3. + (spacing/2.);
                int n_s_bins = int(sigma * 3. / spacing);
                int total_smooth_bins = 1 + (2 * n_s_bins);
                GAUSS<3,float> mygauss(gauss_orig,sigma);
                mygauss.calc_factor_grid(spacing,n_s_bins);

                for (int z=0; z<grs[2]; ++z) {
                    for (int y=0; y<grs[1]; ++y) {
                        for (int x=0; x<grs[0]; ++x) {
                            int gind = x; gind += y*smul1; gind += z*smul2;

                            if (fabs(sg[gind]) < 1.) continue;

                            for (int iz=0; iz<total_smooth_bins; ++iz) {
                                for (int iy=0; iy<total_smooth_bins; ++iy) {
                                    for (int ix=0; ix<total_smooth_bins; ++ix) {
                                        int fgind = ix+x-n_s_bins;
                                        fgind += (iy+y-n_s_bins)*smul1;
                                        fgind += (iz+z-n_s_bins)*smul2;
                                        if (fgind < 0 || fgind >= gs) continue;
                                        fsg[fgind] += mygauss[ix][iy][iz] * sg[gind];
                                    }
                                }
                            }
                        }
                    }
                }
                
                delete[] sg;
                sg = fsg;
            }

            for (int z=0; z<grs[2]; ++z) {
                for (int y=0; y<grs[1]; ++y) {
                    for (int x=0; x<grs[0]; ++x) {
                        int gind = x; gind += y*smul1; gind += z*smul2;
                        if (sg[gind] < min_score) min_score = sg[gind];
                        else if (sg[gind] > max_score) max_score = sg[gind];
                    }
                }
            }

            ofstream f_out1,f_out2,f_out3,f_out4,f_out5;
            string filename1 = lkey + ".acnt";
            string filename2 = lkey + "_good.acnt";
            string filename3 = lkey + "_bad.acnt";
            string filename4 = lkey + "_good_norm.acnt";
            string filename5 = lkey + "_bad_norm.acnt";
            f_out1.open(filename1.c_str());
            if (good_bad) f_out2.open(filename2.c_str());
            if (good_bad) f_out3.open(filename3.c_str());
            if (normalize) f_out4.open(filename4.c_str());
            if (normalize) f_out5.open(filename5.c_str());

            f_out1 << "Atomtype: " << lkey << " Min.-Val.: " << min_score << " Max.-Val.: " << max_score << "\n";
            f_out1 << gmin[ind1] << " " << spacing << " " << gsize[ind1] << "\n";
            f_out1 << gmin[ind2] << " " << spacing << " " << gsize[ind2] << "\n";
            f_out1 << gmin[ind3] << " " << spacing << " " << gsize[ind3] << "\n";
            f_out1 << "90.0 90.0 90.0" << "\n";

            if (good_bad) {
                f_out2 << "Atomtype: " << lkey << " Min.-Val.: " << min_score << " Max.-Val.: " << 0.0 << "\n";
                f_out2 << gmin[ind1] << " " << spacing << " " << gsize[ind1] << "\n";
                f_out2 << gmin[ind2] << " " << spacing << " " << gsize[ind2] << "\n";
                f_out2 << gmin[ind3] << " " << spacing << " " << gsize[ind3] << "\n";
                f_out2 << "90.0 90.0 90.0" << "\n";

                f_out3 << "Atomtype: " << lkey << " Min.-Val.: " << 0.0 << " Max.-Val.: " << max_score << "\n";
                f_out3 << gmin[ind1] << " " << spacing << " " << gsize[ind1] << "\n";
                f_out3 << gmin[ind2] << " " << spacing << " " << gsize[ind2] << "\n";
                f_out3 << gmin[ind3] << " " << spacing << " " << gsize[ind3] << "\n";
                f_out3 << "90.0 90.0 90.0" << "\n";
            }

            if (normalize) {
                f_out4 << "Atomtype: " << lkey << " Min.-Val.: " << 0.0 << " Max.-Val.: " << 1.0 << "\n";
                f_out4 << gmin[ind1] << " " << spacing << " " << gsize[ind1] << "\n";
                f_out4 << gmin[ind2] << " " << spacing << " " << gsize[ind2] << "\n";
                f_out4 << gmin[ind3] << " " << spacing << " " << gsize[ind3] << "\n";
                f_out4 << "90.0 90.0 90.0" << "\n";

                f_out5 << "Atomtype: " << lkey << " Min.-Val.: " << 0.0 << " Max.-Val.: " << 1.0 << "\n";
                f_out5 << gmin[ind1] << " " << spacing << " " << gsize[ind1] << "\n";
                f_out5 << gmin[ind2] << " " << spacing << " " << gsize[ind2] << "\n";
                f_out5 << gmin[ind3] << " " << spacing << " " << gsize[ind3] << "\n";
                f_out5 << "90.0 90.0 90.0" << "\n";
            }

            float tmpval;
            for (int z=0; z<gsize[ind3]; z++) {
                for (int y=0; y<gsize[ind2]; y++) {
                    for (int x=0; x<gsize[ind1]; x++) {
                        int gind = x;
                        if (ind1 == 0) {
                            if (ind2 == 1) { // xyz
                                gind += y*smul1;
                                gind += z*smul2;
                            } else { // xzy
                                gind += z*smul1;
                                gind += y*smul2;
                            }
                        } else if (ind2 == 0) {
                            gind = y;
                            if (ind1 == 1) { // yxz
                                gind += x*smul1;
                                gind += z*smul2;
                            } else { // yzx
                                gind += z*smul1;
                                gind += x*smul2;
                            }
                        } else {
                            gind = z;
                            if (ind1 == 1) { // zxy
                                gind += x*smul1;
                                gind += y*smul2;
                            } else { // zyx
                                gind += y*smul1;
                                gind += x*smul2;
                            }
                        }

                        tmpval = sg[gind];
                        f_out1 << tmpval << "\n";
                        if (tmpval < 0.) {
                            if (good_bad) {
                                f_out2 << tmpval << "\n";
                                f_out3 << 0. << "\n";
                            }
                            if (normalize) {
                                if (min_score < 0.) f_out4 << (tmpval / min_score) << "\n";
                                else f_out4 << 0. << "\n";
                                f_out5 << 0. << "\n";
                            }
                        } else {
                            if (good_bad) {
                                f_out2 << 0. << "\n";
                                f_out3 << tmpval << "\n";
                            }
                            if (normalize) {
                                f_out4 << 0. << "\n";
                                if (max_score > 0.) f_out5 << (tmpval / max_score) << "\n";
                                else f_out5 << 0. << "\n";
                            }
                        }
                    }
                }
            }

            delete[] sg;

            f_out1.close();
            if (good_bad) f_out2.close();
            if (good_bad) f_out3.close();
            if (normalize) f_out4.close();
            if (normalize) f_out5.close();


            pthread_mutex_lock(&wait_mutex);
            ++keys_done;
            if (keys_done == keys_todo) sem_post(&total_keys_worked);
            pthread_mutex_unlock(&wait_mutex);
            if (last_one) {
                sem_wait(&total_keys_worked);
                sem_post(&master_wait);
            }
        }
        pthread_exit(0);
    }
}
using namespace MT;

void remove_dirs(string &filename);

//=================================================================================================
// Deklarationen
//=================================================================================================

class DISTCOUNT {
    private:
        int min_real_bin; //kleinster bin (distanzbasiert)
        int max_real_bin; //groesster bin anhand max_dst
        int min_square_bin;
        int max_square_bin;
        
        float min_dst; //minimale Distanz, die beruecksichtigt wird
        //float max_dst; //maximale Distanz fuer das Scoren
        //float max_square_dst;
        float min_square_dst;
        float target_bin_size;
        float real_bin_size;
        
        //tr1::unordered_map<string,int> key_map;
        
        tr1::unordered_map<string,vector<string> > lig_keys;
        
        GRID<2,float> sgrid; //! sgrid[key_map[key]][bin] = score
        
        //tr1::unordered_map<string,vector<string> > rel_key_map;
        
        tr1::unordered_map<string,float> spec_thresh;
        
        string pot_dir;
        string pdb_file;
        string pro_type;
        string cof_file;
        string wat_file;
        string ref_acnt;
        string acnt_mode;

        //float spacing;
        float shrink;
        
        bool include_metals;
        //bool normalize;
        //bool good_bad;
        
        stl_ptr<STRUCTURE> structure;
        stl_ptr<PARSER> parser;
        stl_ptr<STRUCTURE> cof_structure; //Cofactor
        stl_ptr<PARSER> cof_parser;
        stl_ptr<STRUCTURE> wat_structure; //Wasser
        stl_ptr<PARSER> wat_parser;
        
    public:
        DISTCOUNT(string &pname,string &ptype,string &cname,string &wname,string &pdir,string &racnt,
                  float shrk,bool im,string const& am);
        ~DISTCOUNT();
        
        inline void calc_additional_params();
        inline void load_potentials();
        inline void get_lig_keys();
        inline void load_pdb_file();
        inline void load_map_file(string &map_file);
        inline void generate_rel_key_map();
        void generate();
};


//=================================================================================================
// Definitionen
//=================================================================================================

//!============================================================================================================================
DISTCOUNT::DISTCOUNT(string &pname,string &ptype,string &cname,string &wname,string &pdir,string &racnt,
                     float shrk,bool im,string const& am) {
    pdb_file = pname;
    pro_type = ptype;
    cof_file = cname;
    wat_file = wname;
    pot_dir = pdir;
    ref_acnt = racnt;
    shrink = shrk;
    include_metals = im;
    acnt_mode = am;
    if (am[0] == 'x') MT::ind1 = 0;
    else if (am[0] == 'y') MT::ind1 = 1;
    else if (am[0] == 'z') MT::ind1 = 2;
    if (am[1] == 'x') MT::ind2 = 0;
    else if (am[1] == 'y') MT::ind2 = 1;
    else if (am[1] == 'z') MT::ind2 = 2;
    if (am[2] == 'x') MT::ind3 = 0;
    else if (am[2] == 'y') MT::ind3 = 1;
    else if (am[2] == 'z') MT::ind3 = 2;
    MT::acnt_mode = acnt_mode;
    calc_additional_params();
}

DISTCOUNT::~DISTCOUNT() {
    if (!structure.zero()) structure.kill();
    if (!parser.zero()) parser.kill();
}
//!============================================================================================================================


void DISTCOUNT::calc_additional_params() {
    min_dst = 1.0;
    real_bin_size = 0.01;
    max_dst = 6.0;
    min_square_dst = 1.;
    max_square_dst = 36.;
    target_bin_size = 0.01;
    
    min_real_bin = int(min_dst / real_bin_size);
    max_real_bin = int(max_dst / real_bin_size);
    min_square_bin = int(min_dst * min_dst / target_bin_size);
    max_square_bin = int(max_dst * max_dst / target_bin_size);
}
//!============================================================================================================================


void DISTCOUNT::load_potentials() {
    string ptmp = pot_dir + "/potentials.keys";
    ifstream f_in;
    string row,key_name,last_key;
    f_in.open(ptmp.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << ptmp << "'" << endl;
        structure.kill();
        parser.kill();
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
            structure.kill();
            parser.kill();
            exit(1);
        }
        
        key_map[key_name] = key_num;
        if (key_num > n_diff) n_diff = key_num;
    }
    string shelp = "/potentials_repulsive.bin";
    ptmp = pot_dir + shelp;
    sgrid.load(ptmp.c_str(),true);
    
//    cout << "    -> " << key_map.size() << " different pair potentials loaded from '" << pot_dir << "'\n";
//    cout << "        (" << n_diff+1 << " mapped potentials)" << endl;
    get_lig_keys();
    cout << "     " << lig_keys.size() << " ligand types" << endl;


    gsgrid = new float[sgrid.get_sizes()[0]*sgrid.get_sizes()[1]];
    s0 = sgrid.get_sizes()[0];
    for (int x=0; x<s0; ++x) {
        for (int y=0; y<sgrid.get_sizes()[1]; ++y) {
            int sind = x; sind += y*s0;
            gsgrid[sind] = sgrid[x][y];
        }
    }
}
//!============================================================================================================================


void DISTCOUNT::get_lig_keys() {
    string lk;
    for (tr1::unordered_map<string,int >::iterator it=key_map.begin(); it!=key_map.end(); ++it) {
        string::size_type i1 = it->first.find('_');
        lk = it->first.substr(0,i1);
        if (lig_keys.find(lk) == lig_keys.end()) {
            lig_keys[lk] = vector<string>();
            lig_keys[lk].push_back(lk);
        }
    }
}
//!============================================================================================================================


void DISTCOUNT::load_map_file(string &map_file) {
    ifstream f_in;
    string row,key_name,acnt_name;
    f_in.open(map_file.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << map_file << "'" << endl;
        structure.kill();
        parser.kill();
        exit(1);
    }
    int types_in_file = 0;
    int not_available = 0;
    int mapped_types = 0;
    while (!f_in.eof()) {
        key_name = "@@@";
        acnt_name = "@@@";
        getline(f_in,row);
        istringstream is;
        if (row.size() < 1) continue;
        is.str(row);
        is >> key_name;
        if (key_name[0] == '#') continue;
        if (!is.eof()) is >> acnt_name;
        else acnt_name = key_name;
        if (is.fail() || acnt_name == "@@@") {
            cerr << c_message<cWARNING>("something is wrong with the line '") << row << "' in your map definition file." << endl;
            continue;
        }

        ++types_in_file;
        if (lig_keys.find(key_name) == lig_keys.end()) {
            //cerr << c_message<cWARNING>("no potentials for '") << key_name << "'." << endl;
            ++not_available;
            continue;
        }
        
        if (rel_key_map.find(acnt_name) == rel_key_map.end()) {
            rel_key_map[acnt_name] = vector<string>();
            ++mapped_types;
        }
        for (vector<string>::iterator vkey=lig_keys[key_name].begin(); vkey!=lig_keys[key_name].end(); ++vkey) {
            rel_key_map[acnt_name].push_back(*vkey);
        }
    }
    cout << "     -> " << not_available << " types from your map file are not available" << endl;
    cout << "     -> " << types_in_file - not_available << " types from your map file are mapped to " << mapped_types << " acnt types" << endl;
    f_in.close();
}
//!============================================================================================================================


void DISTCOUNT::generate_rel_key_map() {
    for (tr1::unordered_map<string,vector<string> >::iterator lkey=lig_keys.begin(); lkey!=lig_keys.end(); ++lkey) {
        rel_key_map[lkey->first] = lkey->second;
    }
}
//!============================================================================================================================


void DISTCOUNT::load_pdb_file() {
    string pro_def = pot_dir + "/protein.def";
    structure = new STRUCTURE();
    if (include_metals) parser = new PARSER(&(*structure),1,false,false,true,false);
    else parser = new PARSER(&(*structure),1,false,false,false,false); //!keine Liganden parsen
    
    if (pro_type == "mol") {
        if (!(parser->read_mol2(pdb_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not read ") << pdb_file << endl;
            structure.kill();
            parser.kill();
            exit(1);
        }
        
        cout << "     calculating atom types for '" << pdb_file << "'" << endl;
        
        for (ligands_vec lig=structure->ligands.begin(); lig!=structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,1,true,true);
        }
    } else {
        if (!(parser->read_pdb(pdb_file.c_str()))) {
            cerr << c_message<cERROR>("could not read '") << pdb_file << "' as pdb file" << endl;
            structure.kill();
            parser.kill();
            exit(1);
        }
        
        cout << "     calculating atom types for '" << pdb_file << "'" << endl;
        
        for (chains_vec ct=structure->protein->chains.begin(); ct!=structure->protein->chains.end(); ++ct) {
            //(*ct)->get_atom_typing(1,true,pro_def.c_str(),false,1,true,true); //!Atomtypen der Proteine setzen
            (*ct)->get_atom_typing(1,true,pro_def.c_str(),false,0,true,true);
        }
        
        if (include_metals) {
            for (metals_vec mt=structure->metals.begin(); mt!=structure->metals.end(); ++mt) {
                stl_ptr<LIGAND> nl = new LIGAND();
                stl_ptr<ATOM> na = new ATOM(*(*mt)->atom);
                nl->atoms.push_back(na);
                structure->ligands.push_back(nl);
            }
            for (ligands_vec lig=structure->ligands.begin(); lig!=structure->ligands.end(); ++lig) {
                (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,0,true,true); //!Um die Metalle zu bekommen
            }
        }
    }
    
    structure->name = pdb_file;
    remove_dirs(structure->name);
    
    if (cof_file != "X") {
        cof_structure = new STRUCTURE();
        cof_parser = new PARSER(&(*cof_structure),1,false,true,true,false);
        if (!(cof_parser->read_mol2(cof_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << cof_file << endl;
            exit(1);
        }
        
        cout << "     calculating atom types for '" << cof_file << "'" << endl;
        
        for (ligands_vec lig=cof_structure->ligands.begin(); lig!=cof_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,1,true,true);
        }
    }
    
    if (wat_file != "X") {
        wat_structure = new STRUCTURE();
        wat_parser = new PARSER(&(*wat_structure),1,false,true,true,false);
        if (!(wat_parser->read_mol2(wat_file.c_str(),false,true))) {
            cerr << c_message<cERROR>("could not load ") << wat_file << endl;
            exit(1);
        }
        
        cout << "     calculating atom types for '" << wat_file << "'" << endl;
        
        for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
            (*lig)->get_atom_typing(1,true,pro_def.c_str(),false,1,true,true);
        }
    }
    
}
//!=======================================================================================================


void DISTCOUNT::generate() {
    //!Neue Variante, die schneller laufen sollte:
    //!Ueber alle zu bestimmenden Typen t:
    //!    Ueber alle Proteinatome pa:
    //!        Ueber alle Gitterpunkte gp in einer 12*12*12 A Box um gp:
    //!            pot(t,gp_coord) += pot(t,dst(pa_coord,gp_coord))
    if (pro_type == "mol") {
        for (ligands_vec lig=structure->ligands.begin(); lig!=structure->ligands.end(); ++lig) {
            for (atoms_vec lat=(*lig)->atoms.begin(); lat!=(*lig)->atoms.end(); ++lat) {
                if ((*lat)->sybyl_type == "X") continue;
                grid_atoms.push_back(*lat);
            }
        }
    } else {
        for (chains_vec ct=structure->protein->chains.begin(); ct!=structure->protein->chains.end(); ++ct) {
            for (atoms_vec apro=(*ct)->atoms.begin(); apro!=(*ct)->atoms.end(); ++apro) {
                if ((*apro)->sybyl_type == "X") continue;
                grid_atoms.push_back(*apro);
            }
        }
        if (include_metals) {
            for (ligands_vec lig=structure->ligands.begin(); lig!=structure->ligands.end(); ++lig) {
                for (atoms_vec lat=(*lig)->atoms.begin(); lat!=(*lig)->atoms.end(); ++lat) {
                    if ((*lat)->type == 3) {
                        if ((*lat)->sybyl_type == "X") continue;
                        grid_atoms.push_back(*lat);
                        cout << "     using " << (*lat)->intern_type << " as part of the protein with type "
                        << (*lat)->sybyl_type << endl;
                    }
                }
            }
        }
    }
    if (cof_file != "X") {
        for (ligands_vec lig=cof_structure->ligands.begin(); lig!=cof_structure->ligands.end(); ++lig) {
            for (atoms_vec lat=(*lig)->atoms.begin(); lat!=(*lig)->atoms.end(); ++lat) {
                if ((*lat)->sybyl_type == "X") continue;
                grid_atoms.push_back(*lat);
            }
        }
    }

    //! grid_atoms sollte jetzt alle Proteinatome und zu beruecksichtigenden Metalle enthalten

    if (ref_acnt != "X") {
        //! acnt einlesen und c_min, c_max und spacing bestimmen:
        cout << " --> getting grid parameters from '" << ref_acnt << "'" << endl;
        vec3d<float> c_min; vec3d<float> c_max;
        ifstream f_in;
        string row;
        int n_add;
        f_in.open(ref_acnt.c_str());
        if (!f_in) {
            cerr << c_message<cERROR>("could not read '") << ref_acnt << "'" << endl;
            exit(1);
        }
        getline(f_in,row); row = "";
        getline(f_in,row);
        istringstream is;
        is.str(row);
        is >> c_min[0] >> spacing >> n_add;
        c_max[0] = c_min[0] + (n_add * spacing);
        row = ""; getline(f_in,row);
        is.clear(); is.str(row);
        is >> c_min[1] >> spacing >> n_add;
        c_max[1] = c_min[1] + (n_add * spacing);
        row = ""; getline(f_in,row);
        is.clear(); is.str(row);
        is >> c_min[2] >> spacing >> n_add;
        c_max[2] = c_min[2] + (n_add * spacing);
        f_in.close();
        gmin = c_min;
        gmax = c_max;
    } else {
        gmin = grid_atoms[0]->coord;
        gmax = grid_atoms[0]->coord;
        for (atoms_vec at=grid_atoms.begin(); at!=grid_atoms.end(); ++at) { //zunaechst die Eckpunkte bestimmen
            if ((*at)->coord[0] < gmin[0]) gmin[0] = (*at)->coord[0];
            else if ((*at)->coord[0] > gmax[0]) gmax[0] = (*at)->coord[0];
            if ((*at)->coord[1] < gmin[1]) gmin[1] = (*at)->coord[1];
            else if ((*at)->coord[1] > gmax[1]) gmax[1] = (*at)->coord[1];
            if ((*at)->coord[2] < gmin[2]) gmin[2] = (*at)->coord[2];
            else if ((*at)->coord[2] > gmax[2]) gmax[2] = (*at)->coord[2];
        }

        gmin[0] += shrink; gmax[0] -= shrink;
        gmin[1] += shrink; gmax[1] -= shrink;
        gmin[2] += shrink; gmax[2] -= shrink;
    }

    gmax -= gmin;
    gsize[0] = int((gmax[0] / spacing)+0.5);
    gsize[1] = int((gmax[1] / spacing)+0.5);
    gsize[2] = int((gmax[2] / spacing)+0.5);

    smul1 = gsize[0];
    smul2 = smul1 * gsize[1];

    gs = (1 + gsize[0]) * (1 + gsize[1]) * (1 + gsize[2]);

    acnt_grid = new vec3d<float>[gs];
    
    for (int z=0; z<=gsize[2]; ++z) {
        int gind = z*smul2;
        vec3d<float> hvec(0.,0.,gmin[2] + z*spacing);
        for (int y=0; y<=gsize[1]; ++y) {
            gind += y*smul1;
            hvec[1] = gmin[1] + y*spacing;
            for (int x=0; x<=gsize[0]; ++x) {
                gind += x;
                hvec[0] = gmin[0] + x*spacing;
                acnt_grid[gind] = hvec;
                gind -= x;
            }
            gind -= y*smul1;
        }
    }

    cout << "     using " << spacing << " angstrom grid spacing" << endl;
    cout << "     " << gs << " gridpoints  (" << 1+gsize[0] << "*"
         << 1+gsize[1] << "*" << 1+gsize[2] << ")" << endl;
    cout << "     ";

    //! Erst jetzt noch das Wasser mit dazu nehmen (soll nicht die Gittergroesse mit bestimmen):
    if (wat_file != "X") {
        for (ligands_vec lig=wat_structure->ligands.begin(); lig!=wat_structure->ligands.end(); ++lig) {
            for (atoms_vec lat=(*lig)->atoms.begin(); lat!=(*lig)->atoms.end(); ++lat) {
                if ((*lat)->sybyl_type == "X") continue;
                grid_atoms.push_back(*lat);
            }
        }
    }
    
    //! acnt_grid enthaelt jetzt die Koordinaten aller Gitterpunkte

    // einige Arbeiter erschaffen:
    //cerr << "initializing workers" << endl;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED); // deutlich schneller als JOINABLE !!!
    sem_init(&work_signal,0,0);
    sem_init(&master_wait,0,0);
    sem_init(&total_keys_worked,0,0);
    for (int i=0; i<USE_N_THREADS; ++i) {
        int rc = 0;
        while ((rc = pthread_create(&wt[i],&attr,key_box,NULL))) {}
    }


    for (tr1::unordered_map<string,vector<string> >::iterator lkey=rel_key_map.begin();
                                                              lkey!=rel_key_map.end(); ++lkey) {
        key_stack.push(lkey->first);
    }
    keys_done = 0;
    keys_todo = key_stack.size();
    sem_post(&work_signal);
    sem_wait(&master_wait);



    delete[] acnt_grid;

    cout << endl;
}
//!=======================================================================================================


void remove_dirs(string &filename) {
    string::size_type idx;
    idx = filename.rfind('/');
    if (idx != string::npos) {
        string buf = filename.substr(idx+1);
        filename = buf;
    }
}

void show_header() {
    cout << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n";
    cout << " | 'hotspotsX'     Generating DSX-contour maps (multithreaded version)       |" << "\n";
    cout << " |  author     :   Gerd Neudert                                              |" << "\n";
    cout << " |  supervisor :   Prof. Dr. G. Klebe                  " << c_message<cGREY>("___    _ _") << "            |" << "\n";
    cout << " |  mailto     :   gerd.neudert@gmail.com              " << c_message<cGREY>("))_    )\\`)") << "           |" << "\n";
    cout << " |  version    :   0.61   (17.07.2011)                " << c_message<cGREY>("((_( o ((\\( o") << "          |" << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n" << "\n";
    cout << copyright_string;
    cout << " Call: hotspotsX [-hgmns] -P pro_file [-C cof_file] [-W wat_file]\n";
    cout << "                 [-M maps.def] [-S spacing] [-R shrink] [-D pot_dir]\n"
         << "                 [-A ref.acnt] [-O xyz]\n";
    cout << "=======================================================================" << "\n\n";
}


void show_help() {
    show_header();
    cout << "Type 'hotspotsX -h' to get more help\n" << endl;
}


void show_big_help() {
    show_header();
    cout << " -> pro_file  :  The pdb- or mol2-file you want to create hotspots for.\n";
    cout << "                  From mol2 files everything is taken, from pdb files only\n"
         << "                  ATOM entries will be considered (see cof_file, wat_file\n"
         << "                  and -m option).\n";
    cout << "                  The file should only include those parts of the protein you\n"
         << "                  are interested in. (runtime ~ protein_size^3)\n";
    cout << " -> cof_file  :  Here you can supply a cofactor as a mol2 file. It will be\n"
         << "                  treated as part of the protein.\n";
    cout << " -> wat_file  :  A mol2 file with water molecules (or just O), that will be\n"
         << "                  treated as part of the protein.\n";
    cout << " -> maps.def  :  This file specifies for which atom types you want to generate\n"
         << "                  grids and which maps to combine. Use the first column for\n"
         << "                  the ligand atom types that should be considered and\n"
         << "                  optionally a second column for the acnt-name. If there is\n"
         << "                  the same acnt-name for two (or more) different types, the\n"
         << "                  resulting grid will be the sum of the individual grids.\n";
    cout << "                  Without a maps file you will get an acnt map for each\n"
         << "                  available atomtype of your potential set.\n";
    cout << " -> spacing   :  A float that specifies the grid spacing in Angstrom\n"
         << "                  (default is 0.5)\n";
    cout << " -> shrink    :  The score-grid should be smaller than the full protein-grid\n"
         << "                  to avoid wrong values at the borders. By default the size\n"
         << "                  is 6 angstrom smaller in each direction, but you can change\n"
         << "                  this float value here.\n";
    cout << " -> pot_dir   :  The directory with the potentials to use." << "\n";
    cout << "                  Alternativly you can set the DSX_POTENTIALS environment\n"
         << "                  variable to the pot_dir.\n";
    cout << " -> ref.acnt  :  If you supply a reference acnt-file the grid coordinates will\n"
         << "                  be taken from this file. (So also the spacing etc will be\n"
         << "                  set according to this file!)\n";
    cout << " -> xyz       :  Specify in which order the acnt coordinates should be written.\n"
         << "                  The default is 'xyz', where the x-coords are changing fastest\n"
         << "                  and z-coords slowest. Alternative settings are 'zyx', 'xzy',\n"
         << "                  'zxy', 'yxz' and 'yzx'.\n";
    cout << "                  Please note, that newer versions of PyMOL use 'yxz'!\n";
    cout << " -> -h        :  Prints this help." << "\n";
    cout << " -> -g        :  Write good- and bad-maps (maps with only negative or only\n"
         << "                  positive numbers)\n";
    cout << " -> -m        :  Include metals as part of the protein (for pdb files).\n";
    cout << " -> -n        :  Normalize maps from 0. to 1.\n";
    cout << " -> -s        :  Gaussian smoothing of the final contour maps\n";
    cout << endl;
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
    string cof_file = "X"; //Name des optionalen Kofaktor Files (mol2)
    string wat_file = "X"; //Name des optionalen Wasser Files (mol2)
    string pot_dir = "X"; //Name des Verzeichnis, in dem sie die zu benutzenden Potentiale befinden
    string pro_type = "pdb"; //Dateityp des Protein Files
    string map_file = "X";
    string ref_acnt = "X";
    string acnt_mode = "xyz";
    
    normalize = false;
    bool include_metals = false;
    good_bad = false;
    smooth = false;
    
    spacing = 0.5;
    float shrink = 6.;
    
    //!---------------------------------------------------------------
    //!Kommandozeile parsen:
    for (int i=1; i<argc; ++i) {
        if (s[i][0] == '-') {
            for (int j=1; j<int(s[i].size()); ++j) {
                if (s[i][j] == 'l') {
                    cerr << c_message<cERROR>("the option -l is no longer available since ")
                         << "version 0.50" << endl;
                    cerr << "The default behaviour now is not to write good and bad maps." << endl;
                    cerr << "The flag to switch it on is -g" << endl;
                    show_help();
                    exit(1);
                } else if (s[i][j] == 'm') include_metals = true;
                else if (s[i][j] == 'n') normalize = true;
                else if (s[i][j] == 'g') good_bad = true;
                else if (s[i][j] == 's') smooth = true;
                else if (s[i][j] == 'h') {
                    show_big_help();
                    exit(0);
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
                        map_file = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'S') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -S") << endl;
                        show_help();
                        exit(1);
                    } else if (!string_fu::s2v(s[i],spacing)) {
                        cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                        show_help();
                        exit(1);
                    }
                    break;
                }
                else if (s[i][j] == 'R') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -R") << endl;
                        show_help();
                        exit(1);
                    } else if (!string_fu::s2v(s[i],shrink)) {
                        cerr << c_message<cERROR>(s[i]) << " is not a valid number!" << endl;
                        show_help();
                        exit(1);
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
                else if (s[i][j] == 'A') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -A") << endl;
                        show_help();
                        exit(1);
                    } else {
                        ref_acnt = s[i];
                    }
                    break;
                }
                else if (s[i][j] == 'O') {
                    ++i;
                    if (i == argc) {
                        cerr << c_message<cERROR>("missing argument for -O") << endl;
                        show_help();
                        exit(1);
                    } else {
                        if (s[i] == "xyz" || s[i] == "xzy" || s[i] == "yxz" || 
                            s[i] == "yzx" || s[i] == "zxy" || s[i] == "zyx") acnt_mode = s[i];
                        else {
                            cerr << c_message<cERROR>(s[i]) << " is no valid mode!" << endl;
                            show_help();
                            exit(1);
                        }
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
    
    if (argc < 3) {
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
    if (pro_ext == ".mol2" || pro_ext == ".MOL2") pro_type = "mol";
    else pro_type = "pdb";
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
    
    
    //! 3.) DISTCOUNT-Object erzeugen:
    stl_ptr<DISTCOUNT> mydist = new DISTCOUNT(pro_file,pro_type,cof_file,wat_file,pot_dir,ref_acnt,
                                              shrink,include_metals,acnt_mode);
    
    //! 4.) Potentiale laden:
    cout << " --> loading drugscore potentials" << endl;
    mydist->load_potentials();
    
    //! 5.) PDB laden:
    cout << " --> loading pdb file" << endl;
    mydist->load_pdb_file();
    
    //! 6.) Gucken welche maps gerechnet werden sollen:
    if (map_file != "X") {
        cerr << " --> loading map definition file" << endl;
        mydist->load_map_file(map_file);
    } else {
        cerr << c_message<cWARNING>("you have not specified a maps.def file => hotspots will be calculated for ALL atom types") << endl;
        mydist->generate_rel_key_map();
    }
    
    #if defined (_SPEEDTEST)
    t_total.start_t();
    #endif

    //! 7.) Und jetzt die maps berechnen:
    cout << " --> calculating grids" << endl;
    mydist->generate();

    #if defined (_SPEEDTEST)
    t_total.stop_t();
    #endif

    //! Aufraeumen:
    mydist.kill();
    
    #if defined (_SPEEDTEST)
    cout << "total   = " << t_total.get_time() << endl;
    #endif

    cout << " --> done" << endl;
    if (acnt_mode != "yxz") cout << c_message<cWARNING>("you used acnt_mode '") << acnt_mode << "'. Newer versions of PyMOL use 'yxz' instead (see '-O' option)" << endl;
    
    return 0;
}
