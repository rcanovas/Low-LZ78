/* cdslib - compressed data structures library
    Copyright (C) 2017 Rodrigo Canovas

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

#include <iostream>
#include <unistd.h>
#include "./../include/SLZ78/hlz78.h"
#include "./../include/SLZ78/mhlz78.h"
#include "./../include/SLZ78/ghlz78.h"



using namespace std;

template<class idx_type>
void
decompress(std::string file, std::string out_file) {
    idx_type idx;
    std::ifstream f_in(file, std::ios::in | std::ios::binary);
    if(!f_in) {
        std::cerr << "Failed to open file " << file;
        exit(1);
    }
    idx.load(f_in);
    idx.decompress_file(f_in, out_file);
}

int main(int argc, char* argv[]) {

    if(argc < 2) {
        cout << "Usage: " << argv[0] << " file_name <opt>" << endl;
        cout << "opt: " << endl;
        cout << "-o output_name:  String containing the name of the output file. Default file_name_decom" << endl;
        cout << "-w Index_type. Default = 0" << endl;
        cout << " # | Index_type" << endl;
        cout << "---+--------------------" << endl;
        cout << " 0 | HLZ78 using a map for the displacements" << endl;
        cout << " 1 | HLZ78 using a hash table and a sublayer for displacements" << endl;
        cout << " 2 | MHLZ78 using maps for the displacements" << endl;
        cout << " 3 | MHLZ78 using a hash table and a sublayer for displacements" << endl;
        cout << " 4 | GHLZ78 using maps for the displacements" << endl;
        cout << " 5 | GHLZ78 using a hash table and a sublayer for displacements" << endl;
        return 0;
    }

    string file = argv[1];
    string out_file = file + "_decom";
    uint8_t w = 0;
    double factor = 0.05;
    size_t sigma = 256;

    int c;
    while((c = getopt (argc, argv, "o:w:")) != -1){
        switch (c) {
            case 'o': out_file = optarg;  break;
            case 'w': w = atoi(optarg); break;
            case '?': if(optopt == 'o' || optopt == 'w')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else
                    fprintf(stderr,"Unknown option character `\\x%x'.\n",	optopt);
                return 1;
            default:  abort ();
        }
    }

    //create index
    switch (w) {
        case 0:
            decompress<cdslib::hlz78<> >(file, out_file);
            break;
        case 1:
            decompress<cdslib::hlz78<8, cdslib::hash_D> >(file, out_file);
            break;
        case 2:
            decompress<cdslib::mhlz78<> >(file, out_file);
            break;
        case 3:
            decompress<cdslib::mhlz78<8, cdslib::hash_D> >(file, out_file);
            break;
        case 4:
            decompress<cdslib::ghlz78<> >(file, out_file);
            break;
        case 5:
            decompress<cdslib::ghlz78<8, cdslib::hash_D> >(file, out_file);
            break;
        default:
            cout << "index_type must be a value in [0,5]" << endl;
            break;
    }

    return 0;
}
