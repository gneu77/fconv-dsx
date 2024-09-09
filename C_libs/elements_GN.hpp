
//============================================================================
// elements_GN.hpp -*- C++ -*-; nearly obsolete file, but used by files_GN
//
// Copyright (C) 2006 Gerd Neudert
//
// This file is part of fconv.
// fconv is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------------
//
// author:      Gerd Neudert
// mail:        gerd.neudert@gmail.com
//
// This library is part of the shared libs for fconv, dsx and others.
//
// This library is only used by files_GN for historical reasons
//============================================================================


#ifndef __ELEMENTGN
#define __ELEMENTGN


#include<string>

using namespace std;

//==================================================================================================
//globale Konstanten:
//==================================================================================================

//zum Vergleich mit element-entry:
const string metal_list[] = {"LI" , "BE" , "NA" , "MG" , "AL" , "SI" , " K" , "CA" , "CR" , "MN" , "FE" ,
                             "CO" , "NI" , "CU" , "ZN" , "AG" , "SN" , "CE" , "PT" , "AU" };
//zum Vergleich mit res_name
const string metal_list_ra[] = {" LI" , " BE" , " NA" , " MG" , " AL" , " SI" , "  K" , " CA" , " CR" ,
                                " MN" , " FE" , " CO" , " NI" , " CU" , " ZN" , " AG" , " SN" , " CE" ,
                                " PT" , " AU" };
//zum Vergleich mit name
const string metal_list_la[] = {"LI  " , "BE  " , "NA  " , "MG  " , "AL  " , "SI  " , "K   " , "CA  " ,
                                "CR  " , "MN  " , "FE  " , "CO  " , "NI  " , "CU  " , "ZN  " , "AG  " ,
                                "SN  " , "CE  " , "PT  " , "AU  " };
const int n_metal_list = 20; //Anzahl der explizit bercksichtigten Metalle

const string halogene_list[] = {"F   " , "CL  " , "BR  " , "I   "};
const int n_halogene_list = 4; //Anzahl der bercksichtigten Halogene

//! hier kommen noch vectoren mit zeigern auf halogene, metalle und den rest


//==================================================================================================
//Forward-Deklarationen:
//==================================================================================================

class ELEMENT;


//==================================================================================================
//Deklaration der Klasse element:
//==================================================================================================

class ELEMENT {
    protected:
        string name; //Element-Symbol
        int type; // 2: metall  1: halogen  0: andere
        double weight; // g pro mol
        double cov_radius; //halber Abstand zweier gebundener Atome diesen Typs in Angstrom
        double vdw_radius; //halber Abstand zweier Atome in dichtester Packung (wenn nur vdwWW)
        double pearson_en; //EN nach Pearson in eV
    public:
        ELEMENT();
        ~ELEMENT();
};


#endif //__ELEMENTGN
