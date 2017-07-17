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


#ifndef CDSLIB_MHLZ78_H
#define CDSLIB_MHLZ78_H


#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>
#include <stack>
#include "./../tools.h"
#include "map_D.h"
#include "hash_D.h"

namespace cdslib {

    template<uint8_t t_width = 8, class displacement = map_D>
    class
    mhlz78 {

    public:
        typedef typename sdsl::int_vector<>::size_type                  size_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type    value_type;
        typedef displacement                                            D_type;

    private:
        std::vector<sdsl::int_vector<> *>  H;   //hash tables containing values
        std::vector<D_type *> D;                //structure containing collision information
        std::vector<size_type> M;               //size of each hash table
        std::vector<size_type> P;               //prime number used for each hash table
        size_type number_of_phrases;
        size_type total_size_H;                 //the cumulative sum of the vector M
        size_type sigma;                        // size of the alphabet
        size_type bits_L;


        size_type alpha;                    // constant used for all the tables
        std::vector<size_type> beta;        // mod. mult. inverse of alpha under modulo P[i]
        std::vector<size_type> max_h;       // maximum value of H[i] (indicating that it is empty)


        std::map<value_type, size_type>     alphabet_values;
        sdsl::int_vector<>                  alphabet;


        // auxiliar variables
        uint64_t *buffer;
        uint64_t pos_buffer;
        size_type number_of_collisions;
        size_type max_level;
        size_type cell_used;
        size_type load_limit;
        size_type total_cell_used;
        size_type MAX_N;  //maximum number of phrases
        size_type bits_D;
        double factor;

    public:

        // Empty constructor


        mhlz78() {
            max_level = 0;
            sigma = 0;
            number_of_phrases = 0;
            total_size_H = 0;
        }


        //! Compress file_name into out_file
        //  s : size of the alphabet
        void
        compress_file(std::string file_name, std::string out_file, size_type s,
                      double t_factor = 0.05, size_type bits_d = 0) {
            size_type len_file = 0;
            size_type n = 0;
            bits_D = bits_d;
            factor = t_factor;
            std::ofstream f_out(out_file);
            std::ifstream f_in(file_name, std::ios::in | std::ios::binary);
            if (!f_in) {
                std::cerr << "Failed to open file " << file_name;
                exit(1);
            }
            f_in.seekg(0, f_in.end);
            len_file = f_in.tellg();
            f_in.seekg(0, f_in.beg);
            sigma = s + 2; //+1 to include the symbol of the root and +1 so it is always bigger
            MAX_N = (size_type)std::ceil(len_file / (std::log(1.0 * len_file)/std::log(1.0 * s)));
            bits_L = sdsl::bits::hi((1 + factor) * MAX_N) + 1;
            n = (size_type) std::sqrt(len_file); //could be (1+factor)^2 * sqrt(len_file) 
            M.push_back((size_type) ((1 + factor) * n));
            P.push_back(nearest_prime(M[0] * sigma));
            //std::cout << "len file: " << len_file << "  sigma+2: " << sigma << std::endl;
            //std::cout << "M: " << M[0] << "  P: " << P[0] << std::endl;
            total_size_H = M[0];
            srand(time(NULL)); // initialize random seed*
            alpha = rand() % (P[0] - 1) + 1; //random number in [1, P-1]
            beta.push_back(mod_mul_inv(alpha, P[0]));
            //std::cout << "alpha: " << alpha << "  beta: " << beta[0] << std::endl;
            //need to define values within each H table
            max_h.push_back((P[0] + M[0]) / M[0]);
            //std::cout << "max_h : " << max_h[0] << "  bits_h: " << (sdsl::bits::hi(max_h[0]) + 1) << std::endl;
            H.push_back(new sdsl::int_vector<>(M[0], max_h[0], sdsl::bits::hi(max_h[0]) + 1));
            D.push_back(new D_type(M[0], factor, bits_D));
            load_limit = (size_type)(std::ceil(M[0] * (1.0 / (1.0 + factor))));
            //std::cout << "load limit: " << load_limit << std::endl;
            max_level = 0;
            process_text(f_in, f_out);
            store_data_structure(f_out);
            //std::cout << "total size H: " << total_size_H << std::endl;
        }



        //! Decompress the stored text into out_file. Note that first the data structure must be loaded
        void
        decompress_file(std::istream& in, std::string out_file) {
            in.clear();
            in.seekg(0, in.beg); //make sure that we are at the start of the file
            std::ofstream f_out(out_file);
            size_type h_pos;
            size_type read_size_bits, max_size_bytes, read_space;
            uint64_t end_pos = LoadValue<uint64_t>(in); //read the first variable
            uint64_t current_pos = in.tellg();
            bool need_data = true;
            buffer = cdslib::init_buffer(); //buffer used to write
            pos_buffer = 0;
            uint64_t *reading_buffer;
            uint64_t pos_reading_buffer = 0;
            max_size_bytes = bits_L * (cdslib::BUFFER_SIZE / (8 * bits_L));
            size_type cont = 0;
            while (cont < number_of_phrases) {
                if (need_data) {
                    read_space = end_pos - current_pos;
                    if (read_space > max_size_bytes)
                        read_space = max_size_bytes;
                    read_size_bits = 8 * read_space;
                    reading_buffer = (uint64_t *)(cdslib::LoadValue<char>(in, read_space));
                    current_pos += read_space;
                    need_data = false;
                }
                h_pos = cdslib::get_field(reading_buffer, pos_reading_buffer, pos_reading_buffer + bits_L - 1);
                extract(h_pos, f_out);
                pos_reading_buffer += bits_L;
                if (pos_reading_buffer == read_size_bits) {
                    delete [] reading_buffer;
                    pos_reading_buffer = 0;
                    need_data = true;
                }
                ++cont;
            }
            if (pos_buffer != 0)
                cdslib::SaveValue(f_out, (char *)buffer, (pos_buffer + 7) / 8);
            delete [] buffer;
            if (pos_reading_buffer > 0)
                delete [] reading_buffer;
        }

        //! Loads the data structure
        void
        load(std::istream& in) {
            uint64_t start_position = LoadValue<uint64_t>(in);
            in.seekg(start_position);
            sdsl::read_member(number_of_phrases, in);
            sdsl::read_member(sigma, in);
            sdsl::read_member(bits_L, in);
            sdsl::read_member(alpha, in);
            sdsl::read_member(total_size_H, in);
            sdsl::read_member(max_level, in);

            M = std::vector<size_type>(max_level + 1);
            P = std::vector<size_type>(max_level + 1);
            beta = std::vector<size_type>(max_level + 1);
            max_h = std::vector<size_type>(max_level + 1);
            H = std::vector<sdsl::int_vector<> *>(max_level + 1);
            D = std::vector<D_type *>(max_level + 1);
            for (size_type i = 0; i <= max_level; ++i) {
                sdsl::read_member(M[i], in);
                sdsl::read_member(P[i], in);
                sdsl::read_member(beta[i], in);
                sdsl::read_member(max_h[i], in);
                {
                    H[i] = new sdsl::int_vector<>();
                    sdsl::int_vector<> B;
                    B.load(in); // load B
                    load_vector_selection((*H[i]), B, max_h[i], in);
                }
                D[i] = new D_type();
                D[i]->load(in);
            }
            alphabet.load(in); //Load map_alphabet
        }


    private:

        void
        store_data_structure(std::ostream& out) {
            uint64_t start_position = out.tellp(); //get actual position in the file
            std::cout << "Finish storing L at:" << start_position << std::endl;
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(nullptr, "", sdsl::util::class_name(*this));
            size_type size_of_struc = 0, size_h = 0, size_d = 0;
            size_of_struc += write_member(number_of_phrases, out, child, "num phrases");
            size_of_struc += write_member(sigma, out, child, "sigma");
            size_of_struc += write_member(bits_L, out, child, "bits L");
            size_of_struc += write_member(alpha, out, child, "alpha");
            size_of_struc += write_member(total_size_H, out, child, "cumulative size of the table");
            size_of_struc += write_member(max_level, out, child, "max_level");
            for (size_type i = 0; i <= max_level; ++i) {
                size_of_struc += write_member(M[i], out, child, "size table");
                size_of_struc += write_member(P[i], out, child, "Prime number");
                size_of_struc += write_member(beta[i], out, child, "beta value");
                size_of_struc += write_member(max_h[i], out, child, "max_h value");

                sdsl::int_vector<> B(M[i], 0, 1);
                size_type count = 0;
                for(size_type j = 0; j < M[i]; ++j) {
                    if ((*H[i])[j] != max_h[i]) {
                        B[j] = 1;
                        ++count;
                    }
                }
                size_of_struc += B.serialize(out, child, "bitmap B"); //Save B
                size_h += write_vector_selection((*H[i]), B, count, out);
                size_d += (*D[i]).serialize(out, child, "D structure");
            }
            size_of_struc += size_d + size_h;
            //std::cout << "total length H: " << total_size_H << std::endl;
            //std::cout << "H: " << size_h << " bytes" << std::endl;
            //std::cout << "D: " << size_d << " bytes" << std::endl;
            alphabet.serialize(out, child, "map alphabet"); //Save map_alphabet
            out.seekp(0,std::ios::beg);
            SaveValue(out, start_position);
        }


        void
        process_text(std::istream& in, std::ostream& out) {
            SaveValue(out, (int64_t)0); //save space to store the start of the data structure
            size_type current_node_pos, node_value;
            size_type count_sigma = 0;
            size_type current_index = 0;
            number_of_phrases = number_of_collisions = 0;
            cell_used = 0;
            total_cell_used = 0;
            alphabet = sdsl::int_vector<>(sigma, 0, t_width);
            pos_buffer = 0;
            buffer = cdslib::init_buffer();
            (*H[current_index])[0] = 0; //insert the root to the hash at position 0
            ++cell_used;
            current_node_pos = 0;
            char c;
            while (in.get(c)) {
                if (count_sigma < sigma) {
                    auto it = alphabet_values.find(c);
                    if (it == alphabet_values.end()) {
                        alphabet_values[c] = count_sigma + 1;
                        alphabet[count_sigma] = (size_type)c;
                        ++count_sigma;
                    }
                }
                node_value = current_node_pos * sigma + alphabet_values[c];
                current_node_pos = insert(node_value, out, current_index);
            }
            if(current_node_pos != 0) { //last phrase
                ++number_of_phrases;
                cdslib::set_field(buffer, pos_buffer, pos_buffer + bits_L, current_node_pos);
                pos_buffer += bits_L;
            }
            //std::cout << "number of levels: " << M.size() << std::endl;
            //std::cout << "alphabet_size: " << alphabet_values.size() << std::endl;
            //std::cout << "number of phrases: " << number_of_phrases << std::endl;
            //std::cout << "number of collisions: " << number_of_collisions << std::endl;
            //std::cout << "number of levels: " << max_level + 1 << std::endl;
            if (pos_buffer != 0)
                cdslib::SaveValue(out, buffer, (pos_buffer + 63) / 64);
            delete [] buffer;
        }


        //! Insert the h_value at h_pos if possible. Otherwise check using linear
        // probing if the value already exist or need to be inserted. Returns
        // the position of the node where the next search should start.
        size_type
        insert(size_type node_value, std::ostream& out, size_type &current_index) {
            size_type d = 0;
            bool not_found = false;
            size_type h_value = cdslib::my_mod2(alpha, node_value, P[current_index]);
            size_type h_pos = h_value % M[current_index];
            size_type original_h_pos = h_pos;
            h_value = h_value / M[current_index];
            if ((*H[current_index])[h_pos] == max_h[current_index]) { //is empty
                not_found = true;
            }
            else {
                while ((*H[current_index])[h_pos] != max_h[current_index]) {
                    if (verify(original_h_pos, h_pos, h_value, current_index)) {//if already exist
                        if (current_index != 0)
                            h_pos += (1 << (current_index - 1)) * M[0];
                        return h_pos;
                    }
                    ++d;
                    if (h_pos == M[current_index] - 1)  h_pos = 0;
                    else                                h_pos = h_pos + 1;
                }
                not_found = true;
            }
            if (not_found) {
                if (current_index == max_level) {
                    if (cell_used < load_limit) {
                        (*H[current_index])[h_pos] = h_value;
                        D[current_index]->setValue(h_pos, d);
                        ++number_of_phrases;
                        if (d > 0)
                            ++ number_of_collisions;
                        if (current_index != 0)
                            h_pos += (1 << (current_index - 1)) * M[0];
                        cdslib::set_field(buffer, pos_buffer, pos_buffer + bits_L, h_pos);
                        pos_buffer += bits_L;
                        cdslib::check_buffer(buffer, pos_buffer, out);
                        ++cell_used;
                        current_index = 0;
                        return 0; //come back to the root
                    }
                    else {
                        create_new_level();
                        current_index = current_index + 1;
                        return insert(node_value, out, current_index);
                    }
                }
                else {
                    current_index = current_index + 1;
                    return insert(node_value, out, current_index);
                }
            }
        }

        void
        create_new_level() {
            size_type level_size = (1 << max_level) * M[0]; // 2^(max_level)*M[0]
            total_cell_used += cell_used;
            size_type max_left_phrases = (size_type)std::ceil((MAX_N - total_cell_used) * (1.0 + factor));
            level_size = (size_type)std::min(level_size, max_left_phrases);
            ++max_level;
            M.push_back(level_size);
            total_size_H += level_size;
            P.push_back(nearest_prime(total_size_H * sigma));
            beta.push_back(mod_mul_inv(alpha, P[max_level]));
            max_h.push_back((P[max_level] + M[max_level]) / M[max_level]); //need to define values within the H table
            H.push_back(new sdsl::int_vector<>(M[max_level], max_h[max_level], sdsl::bits::hi(max_h[max_level]) + 1));
            D.push_back(new D_type(M[max_level], factor, bits_D));
            load_limit = (size_type)(std::ceil(M[max_level] * (1.0 / (1.0 + factor))));
            cell_used = 0;
        }

        //! Returns true if the value stored at pos in the Hash contains the same
        // information than h_value
        bool
        verify(size_type original_h_pos, size_type h_pos, size_type h_value, size_type &current_index) {
            if (h_pos == 0 and current_index == 0)  //the root
                return false;
            size_type pos = h_pos; //position in the hash table
            size_type mod_P_value = 0, n_value = 0;
            if ((*H[current_index])[h_pos] != h_value)
                return false;
            //update position in case depending of the collisions
            pos = (pos + M[current_index] - D[current_index]->getValue(h_pos)) % M[current_index]; //update position in case depending of the collisions
            if (pos != original_h_pos)
                return false;
            return true; //pos == original_h_pos  and H[h_pos] == h_value
        }


        //!Extracts the value path from h_pos to the root into the out file
        void
        extract(size_type h_pos, std::ostream& out) {
            std::stack<value_type > phrase;
            size_type local_h_pos, current_level;
            size_type current_h_value, actual_h_pos;
            size_type mod_P_value, node_value;
            while (h_pos != 0) {
                if (h_pos < M[0]) {
                    current_level = 0;
                    local_h_pos = h_pos;
                }
                else {
                    current_level = sdsl::bits::hi((size_type)(h_pos/M[0])) + 1;
                    local_h_pos = h_pos - (1 << (current_level - 1)) * M[0];
                }
                current_h_value = (*H[current_level])[local_h_pos];
                actual_h_pos = local_h_pos + M[current_level]
                               - D[current_level]->getValue(local_h_pos);
                actual_h_pos = actual_h_pos % M[current_level];
                mod_P_value = current_h_value * M[current_level] + actual_h_pos;
                node_value = cdslib::my_mod2(mod_P_value, beta[current_level], P[current_level]); //(mod_P_value * beta) % P;
                h_pos = node_value / sigma;
                phrase.push((value_type)(node_value - h_pos * sigma));
            }
            while (!phrase.empty()) { //write to the buffer
                cdslib::set_field(buffer, pos_buffer, pos_buffer + t_width, alphabet[phrase.top() - 1]);
                phrase.pop();
                pos_buffer += t_width;
                cdslib::check_buffer(buffer, pos_buffer, out);
            }
        }

    };

}

#endif //CDSLIB_HLZ78_H
