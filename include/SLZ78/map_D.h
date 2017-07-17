/* Low-LZ78
   Copyright (C) 2016-2017 Rodrigo Canovas
    
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


#ifndef CDSLIB_MAP_D_H
#define CDSLIB_MAP_D_H

#include "./../tools.h"

namespace cdslib {

    class map_D {
    public:
        typedef typename sdsl::int_vector<>::size_type  size_type;

    private:
        sdsl::int_vector<> D;               //array containing collision information
        std::map<size_type, size_type> E;   //special case collisions
        size_type max_d;                    // maximum value stored in D

    public:

        // Empty constructor
        map_D() { }

        // (t_d_bits = 3 bits original code)
        map_D(size_type M, double factor, size_type t_d_bits = 0) {
            size_type bits_D = t_d_bits;
             if (bits_D == 0) {
                 max_d = (size_type) (1.0 / factor);
                 bits_D = sdsl::bits::hi(max_d) + 1;
             }
             else
                 max_d = (1 << bits_D) - 1;
            D = sdsl::int_vector<>(M, 0, bits_D);
        }

        void
        swap(map_D& md) {
            if (this != &md) {
                D.swap(md.D);
                std::swap(max_d, md.max_d);
                E.swap(md.E);
            }
        }


        size_type
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0, e_bytes = 0;
            written_bytes += write_member(max_d, out, child, "maximum value in D");
            written_bytes += D.serialize(out, child, "D array");
            {//Save E
                size_type length_E = E.size(), cont = 0;
                //std::cout << "number of special cases: " << length_E << std::endl;
                sdsl::int_vector<> elements(length_E, 0, 8*sizeof(size_type));
                for(auto it = E.begin(); it != E.end(); it++) {
                    elements[cont] = it->first;
                    cont ++;
                }
                e_bytes = elements.serialize(out, child, "map keys");
                cont = 0;
                for(auto it = E.begin(); it != E.end(); it++) {
                    elements[cont] = it->second;
                    cont ++;
                }
                e_bytes += elements.serialize(out, child, "map values");
            }
            //std::cout << "Space used by special cases: " << e_bytes << " bytes" << std::endl;
            written_bytes += e_bytes;
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void
        load(std::istream& in) {
            sdsl::read_member(max_d, in);
            D.load(in);
            sdsl::int_vector<> keys;
            sdsl::int_vector<> values;
            keys.load(in);
            values.load(in);
            for(size_type i = 0; i < keys.size(); ++i)
                E[keys[i]] = values[i];
        }


        size_type
        getValue(size_type pos) {
            if (D[pos] == max_d) {
                auto it = E.find(pos); // assume that always find it
                return it->second;
            }
            else
                return D[pos];
        }

        void
        setValue(size_type pos, size_type d) {
            if (d >= max_d) {
                D[pos] = max_d;
                E[pos] = d;
            }
            else
                D[pos] = d;
        }

    };


}

#endif //CDSLIB_MAP_D_H
