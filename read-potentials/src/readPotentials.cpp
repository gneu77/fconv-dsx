#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<map>
#include<tr1/unordered_map>
#include"grid_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"
#include"string_fu_GN.hpp"

using namespace std;
using namespace TXT;

const string version_string = "0.10";
const string version_date = "09.09.2024";
const string copyright_string =
" Copyright (C) 2024 Gerd Neudert\n\
  Usage is free without any limitations.\n\
  There is NO warranty, not even for MERCHANTABILITY or FITNESS\n\
  FOR A PARTICULAR PURPOSE.\n\n";

//!============================================================================================================================
void show_header() {
    cout << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n";
    cout << " | 'readPotentials' : Parse binary DSX potentials and output readable format |" << "\n";
    cout << " |  author          : Gerd Neudert                                           |" << "\n";
    cout << " |  version         : " << version_string << "   (" << version_date << ")                                    |" << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n" << "\n";
    cout << copyright_string;
    cout << " => Call:" << "\n";
    cout << "  readPotentials [-h] -D pot_dir\n";
    cout << "========================================================================" << "\n" << "\n";
}

void show_help() {
    show_header();
    cout << "Type 'readPotentials -h' to get more help" << "\n" << endl;
}

void show_big_help() {
    show_header();
    cout << "  pot_dir     :  The directory with the potentials to use.\n";
    cout << "\n Output will be on the console";
    cout << endl;
}

int main(int argc,char *argv[]) {
    // Spaghetti code taken from drugscoreX.cpp:

    vector<string> s;
    s.resize(argc);
    for (int i=1; i<argc; ++i) {
        s[i].assign(argv[i]);
    }
    for (int i=1; i<argc; ++i) { //Leerzeichen entfernen
        ostringstream os;
        istringstream is;
        is.clear(); is.str(s[i]);
        string tmp; is >> tmp;
        os << tmp; s[i] = os.str();
    }
    
    string pot_dir = "X";
    
    for (int i=1; i<argc; ++i) {
        if (s[i][0] == '-') {
            for (int j=1; j<int(s[i].size()); ++j) {
                if (s[i][j] == 'h') {
                    show_big_help();
                    exit(0);
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
    
    if (pot_dir == "X") {
        cerr << c_message<cERROR>("you have to specify the directory with the set of potentials to use") << endl;
        show_help();
        exit(1);
    }
    
    // read potentials:
    string ptmp = pot_dir + "/potentials.keys";
    ifstream f_in;
    string row,key_name,last_key;
    f_in.open(ptmp.c_str());
    if (!f_in) {
        cerr << c_message<cERROR>("could not read '") << ptmp << "'" << endl;
        exit(1);
    }
    getline(f_in,row);
    f_in.close();
    istringstream is;
    is.str(row);
    int key_num;
    last_key = "*";
    int n_diff = 0;
    tr1::unordered_map<string,int> key_map;
    tr1::unordered_map<int,vector<string>> key_map_inv;
    while (!is.eof()) {
        key_name = last_key;
        is >> key_name;
        if (key_name == last_key) break;
        key_num = -1;
        is >> key_num;
        if (key_num == -1) {
            cerr << c_message<cERROR>(ptmp) << " is not valid (only potentials newer than 08/2009 are valid)" << endl;
            exit(1);
        }
        key_map[key_name] = key_num;
        vector<string> tv;
        string_fu::mysplit(key_name,tv,'_');
        key_map_inv[key_num] = tv;
        if (key_num > n_diff) n_diff = key_num;
    }
    cout << "<potential keys> " << ptmp << "\n";

    string shelp = "/potentials_repulsive.bin";
    ptmp = pot_dir + shelp;
    GRID<2,float> sgrid; //! sgrid[key_map[key]][bin] = score  // the pair-potentials
    sgrid.load(ptmp.c_str(),true);
    cout << "<repulsive pair potentials> " << ptmp << "\n";

    // output:
    int nKeys = sgrid.get_sizes()[0];
    int nBins = sgrid.get_sizes()[1];
    cout << "<output format> protein-type | ligand-type | distance | potential-value" << "\n";
    for (int keyIndex=0; keyIndex<nKeys; ++keyIndex) {
        for (int bin=0; bin<nBins; ++bin) {
            cout << key_map_inv[keyIndex][0] << " | " << key_map_inv[keyIndex][1] << " | " << (bin/100.) << " | " << sgrid[keyIndex][bin] << "\n";
        }
    }
    cout << endl;
    return 0;
}
