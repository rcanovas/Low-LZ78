/* Low-LZ78 - compressed data structures library
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
compress_text(std::string file, std::string out_file, double factor, size_t sigma, size_t d_bits) {
    auto idx = idx_type();
    idx.compress_file(file, out_file, sigma, factor, d_bits);
		//std::cout << "Size Index: " << sdsl::size_in_bytes(idx) << " bytes" << std::endl;
}

int main(int argc, char* argv[]) {

    if(argc < 2) {
        cout << "Usage: " << argv[0] << " file_name <opt>" << endl;
        cout << "opt: " << endl;
        cout << "-o output_name:  String containing the name of the output file. Default file_name.hlz78" << endl;
        cout << "-s sigma: Size of the alphabet of the input text. Default s = 256" << endl;
        cout << "-f factor: Integer value indicating the overload factor used  for the Hash table. Default factor = 5%" << endl;
        cout << "-d D_bits: Number of bits used to store each collision value. Default d=0 and compute it internally" << endl;
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
    string out_file = file;
    uint8_t w = 0;
    double factor = 0.05;
    size_t sigma = 256;
    size_t d_bits = 0;

    int c;
    while((c = getopt (argc, argv, "o:w:f:s:d:")) != -1){
        switch (c) {
            case 'o': out_file = optarg;  break;
            case 'w': w = atoi(optarg); break;
            case 's': sigma = atoi(optarg); break;
            case 'f': factor = 1.0 * atoi(optarg) / 100.0; break;
            case 'd': d_bits = atoi(optarg); break;
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
            out_file += ".hlz78";
            compress_text<cdslib::hlz78<> >(file, out_file, factor, sigma, d_bits);
            break;
        case 1:
            out_file += ".hlz78_hash";
            compress_text<cdslib::hlz78<8, cdslib::hash_D> >(file, out_file, factor, sigma, d_bits);
            break;
        case 2:
            out_file += ".mhlz78";
            compress_text<cdslib::mhlz78<> >(file, out_file, factor, sigma, d_bits);
            break;
        case 3:
            out_file += ".mhlz78_hash";
            compress_text<cdslib::mhlz78<8, cdslib::hash_D> >(file, out_file, factor, sigma, d_bits);
            break;
        case 4:
            out_file += ".ghlz78";
            compress_text<cdslib::ghlz78<> >(file, out_file, factor, sigma, d_bits);
            break;
        case 5:
            out_file += ".ghlz78_hash";
            compress_text<cdslib::ghlz78<8, cdslib::hash_D> >(file, out_file, factor, sigma, d_bits);
            break;
        default:
            cout << "index_type must be a value in [0,5]" << endl;
            break;
    }

    return 0;
}
