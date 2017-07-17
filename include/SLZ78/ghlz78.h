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


#ifndef CDSLIB_GMHLZ78_H
#define CDSLIB_GMHLZ78_H

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
    ghlz78 {

    public:
        typedef typename sdsl::int_vector<>::size_type size_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type value_type;
        typedef displacement D_type;

    private:
        sdsl::int_vector<> H;               //hash table containing values
        D_type D;                           //structure containing collision information
        size_type number_of_phrases;
        size_type M;                        //size of the hash table H
        size_type P;                        // prime number
        size_type alpha;                    // constant used
        size_type beta;                     // mod. mult. inverse of alpha under modulo P
        size_type sigma;                    // size of the alphabet
        size_type max_h;                    // maximum value of H (indicating that it is empty)
        size_type bits_L;


        size_type cell_used;
        size_type load_limit;
        size_type bits_d;
        size_type MAX_M;
        double factor;

        std::map<value_type, size_type> alphabet_values;
        sdsl::int_vector<> alphabet;

        uint64_t *buffer;
        uint64_t pos_buffer;
        std::string out_file_name;
        size_type number_of_collisions;
        size_type max_length;

    public:

        // Empty constructor
        ghlz78() {
            M = P = alpha = beta = sigma = 0;
            number_of_phrases = max_h = max_length = 0;
            load_limit = cell_used = 0;
        }

        void
        swap(ghlz78& gh) {
            if (this != &gh) {
                H.swap(gh.H);
                D.swap(gh.D);
                std::swap(M, gh.M);
                std::swap(MAX_M, gh.MAX_M);
                std::swap(sigma, gh.sigma);
                std::swap(alpha, gh.alpha);
                std::swap(bits_d, gh.bits_d);
                std::swap(factor, gh.factor);
                std::swap(P, gh.P);
                std::swap(beta, gh.beta);
                std::swap(max_h, gh.max_h);
                std::swap(number_of_phrases, gh.number_of_phrases);
                std::swap(number_of_collisions, gh.number_of_collisions);
                std::swap(load_limit, gh.load_limit);
                std::swap(cell_used, gh.cell_used);
                std::swap(bits_L, gh.bits_L);
                std::swap(buffer, gh.buffer);
                std::swap(pos_buffer, gh.pos_buffer);
            }
        }

        //! Compress file_name into out_file
        //  s : size of the alphabet
        void
        compress_file(std::string file_name, std::string out_file, size_type s,
                      double _factor = 0.05, size_type _bits_d = 0) {
            size_type len_file = 0;
            size_type n = 0, bits_H = 0;
            out_file_name = out_file;
            std::ofstream f_out(out_file);
            std::ifstream f_in(file_name, std::ios::in | std::ios::binary);
            if (!f_in) {
                std::cerr << "Failed to open file " << file_name;
                exit(1);
            }
            f_in.seekg(0, f_in.end);
            len_file = f_in.tellg();
            f_in.seekg(0, f_in.beg);
            //initialize H and D, and variables P,M,alpha,beta, and sigma
            sigma = s + 2; //+1 to include the symbol of the root and +1 so it is always bigger
            factor = _factor;
            n = (size_type) std::ceil(len_file / (std::log(1.0 * len_file) / std::log(1.0 * s))); //worse case
            MAX_M = (1 + factor) * n;
            bits_L = sdsl::bits::hi(MAX_M) + 1;
            n = (size_type) std::sqrt(len_file); //best case-->initial case   (could be (1+factor)^2 * sqrt(len_file))
            M = (size_type) ((1 + factor) * n);
            P = nearest_prime(M * sigma);
            //std::cout << "len file: " << len_file << "  sigma+2: " << sigma << std::endl;
            //std::cout << "M: " << M << "  P: " << P << std::endl;
            srand(time(NULL)); // initialize random seed
            alpha = rand() % (P - 1) + 1; //random number in [1, P-1]
            beta = mod_mul_inv(alpha, P);
            //std::cout << "alpha: " << alpha << "  beta: " << beta << std::endl;
            max_h = (P + M) / M;
            bits_H = sdsl::bits::hi(max_h) + 1; //aprox log_2(sigma) in the original code log_2(sigma+2)
            H = sdsl::int_vector<>(M, max_h, bits_H); //(size, default value, bits_per_symbol)
            bits_d = _bits_d;
            D = D_type(M, factor, bits_d);
            load_limit = (size_type)(std::ceil(M * (1.0 / (1.0 + factor))));
            //maximum number of nodes possible (worse case)
            process_text(f_in, f_out);
            store_data_structure(f_out);
        }

        //!Extracts the value path from h_pos to the root into the out file
        std::vector<value_type>
        extract(size_type h_pos) {
            std::vector<value_type> phrase;
            size_type current_h_value, actual_h_pos, old_h_pos;
            size_type mod_P_value, node_value, d_value;
            while (h_pos != 0) {
                old_h_pos = h_pos; //this is new
                current_h_value = H[h_pos];
                actual_h_pos = h_pos;
                d_value = D.getValue(h_pos);
                if (h_pos >= d_value)
                    actual_h_pos -= d_value;
                else
                    actual_h_pos += M - d_value;
                mod_P_value = current_h_value * M + actual_h_pos;
                node_value = cdslib::my_mod2(mod_P_value, beta, P);//(mod_P_value * beta) % P;
                h_pos = node_value / sigma;
                phrase.push_back((value_type) (node_value - h_pos * sigma));
            }
            return phrase;
        }

         //! Insert the h_value at h_pos if possible. Otherwise check using linear
        // probing if the value already exist or need to be inserted. Returns
        // the position of the node where the next search should start.
        size_type
        insert(size_type node_value, std::ostream& out) {
             size_type h_value = cdslib::my_mod2(alpha, node_value, P);
             size_type h_pos = h_value % M;
             h_value = h_value / M;
             size_type d = 0;
             size_type original_h_pos = h_pos;
             bool not_found = false;
             if (H[h_pos] == max_h)//is empty so need to insert
                 not_found = true;
             else {
                 while (H[h_pos] != max_h) {
                     if (verify(original_h_pos, h_pos, h_value))//if already exist
                         return h_pos;
                     ++d;
                     if (h_pos == M - 1) h_pos = 0;
                     else h_pos = h_pos + 1;
                 }
                 not_found = true;
                 ++number_of_collisions;
             }
             if (not_found) {
                 H[h_pos] = h_value;
                 ++number_of_phrases;
                 D.setValue(h_pos, d);
                 cdslib::set_field(buffer, pos_buffer, pos_buffer + bits_L, h_pos);
                 pos_buffer += bits_L;
                 cdslib::check_buffer(buffer, pos_buffer, out);
                 ++cell_used;
                 if (cell_used == load_limit)
                     grow_tables(out);
                 return 0; //come back to the root
             }
         }


       //! Decompress the stored text into out_file. Note that first the data structure must be loaded
        void
        decompress_file(std::istream &in, std::string out_file) {
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
           max_length = 0;
           while (cont < number_of_phrases) {
               if (need_data) {
                   read_space = end_pos - current_pos;
                   if (read_space > max_size_bytes)
                       read_space = max_size_bytes;
                   read_size_bits = 8 * read_space;
                   reading_buffer = (uint64_t *) (cdslib::LoadValue<char>(in, read_space));
                   current_pos += read_space;
                   need_data = false;
               }
               h_pos = cdslib::get_field(reading_buffer, pos_reading_buffer, pos_reading_buffer + bits_L - 1);
               extract(h_pos, f_out);
               pos_reading_buffer += bits_L;
               if (pos_reading_buffer == read_size_bits) {
                   delete[] reading_buffer;
                   pos_reading_buffer = 0;
                   need_data = true;
               }
               ++cont;
           }
           if (pos_buffer != 0)
               cdslib::SaveValue(f_out, (char *) buffer, (pos_buffer + 7) / 8);
           delete[] buffer;
           if (pos_reading_buffer > 0)
               delete[] reading_buffer;
       }

        //! Loads the data structure
        void
        load(std::istream &in) {
            uint64_t start_position = LoadValue<uint64_t>(in);
            in.seekg(start_position);
            size_type *parameters;
            parameters = LoadValue<size_type>(in, 8);
            number_of_phrases = parameters[0];
            M = parameters[1];
            P = parameters[2];
            alpha = parameters[3];
            beta = parameters[4];
            sigma = parameters[5];
            max_h = parameters[6];
            bits_L = parameters[7];
            sdsl::int_vector<> B;
            B.load(in); // load B
            load_vector_selection(H, B, max_h, in);
            D.load(in);//Load D
            alphabet.load(in); //Load map_alphabet
            delete[] parameters;
        }

    private:

        ghlz78(size_type _MAX_M, size_type _M, size_type _sigma,
               size_type _alpha, size_type _bits_d, double _factor, size_type _bits_L) {
            M = _M;
            MAX_M = _MAX_M;
            sigma = _sigma;
            alpha = _alpha;
            bits_d = _bits_d;
            factor = _factor;
            bits_L = _bits_L;
            P = nearest_prime(M * sigma);
            beta = mod_mul_inv(alpha, P);
            max_h = (P + M) / M;
            H = sdsl::int_vector<>(M, max_h, sdsl::bits::hi(max_h) + 1);
            H[0] = 0; //insert the root to the hash at position 0
            cell_used = 1;
            D = D_type(M, factor, bits_d);
            load_limit = (size_type)(std::ceil(M * (1.0 / (1.0 + _factor))));
            number_of_phrases = number_of_collisions = 0;
            pos_buffer = 0;
            buffer = cdslib::init_buffer();
        }


        void
        store_data_structure(std::ostream &out) {
            uint64_t start_position = out.tellp(); //get actual position in the file
            //std::cout << "Finish storing L at:" << start_position << std::endl;
            size_type size_of_struc = 0;
            size_type parameters[8];
            parameters[0] = number_of_phrases;
            parameters[1] = M;
            parameters[2] = P;
            parameters[3] = alpha;
            parameters[4] = beta;
            parameters[5] = sigma;
            parameters[6] = max_h;
            parameters[7] = bits_L;
            SaveValue(out, parameters, 8);
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(nullptr, "",
                                                                               sdsl::util::class_name(*this));
            sdsl::int_vector<> B(M, 0, 1);
            size_type count = 0;
            for(size_type i = 0; i < M; ++i) {
                if (H[i] !=  max_h) {
                    B[i] = 1;
                    ++count;
                }
            }
            size_of_struc = B.serialize(out, child, "bitmap B"); //Save B
            write_vector_selection(H, B, count, out);

            size_of_struc = D.serialize(out, child, "D Array"); //Save D
            alphabet.serialize(out, child, "map alphabet"); //Save map_alphabet
            out.seekp(0, std::ios::beg);
            SaveValue(out, start_position);
        }

        void
        process_text(std::istream &in, std::ostream &out) {
            SaveValue(out, (int64_t) 0); //save space to store the start of the data structure
            size_type current_node_pos, h_pos, h_value, node_value;
            size_type count_sigma = 0;
            number_of_phrases = number_of_collisions = 0;
            alphabet = sdsl::int_vector<>(sigma, 0, t_width);
            pos_buffer = 0;
            buffer = cdslib::init_buffer();
            H[0] = 0; //insert the root to the hash at position 0
            cell_used = 1;
            current_node_pos = 0;
            char c;
            while (in.get(c)) {
                if (count_sigma < sigma) {
                    auto it = alphabet_values.find(c);
                    if (it == alphabet_values.end()) {
                        alphabet_values[c] = count_sigma + 1;
                        alphabet[count_sigma] = (size_type) c;
                        ++count_sigma;
                    }
                }
                node_value = current_node_pos * sigma + alphabet_values[c];
                current_node_pos = insert(node_value, out);
            }
            if (current_node_pos != 0) { //last phrase
                ++number_of_phrases;
                cdslib::set_field(buffer, pos_buffer, pos_buffer + bits_L, current_node_pos);
                pos_buffer += bits_L;
            }
            //std::cout << "alphabet_size: " << alphabet_values.size() << std::endl;
            //std::cout << "number of phrases: " << number_of_phrases << std::endl;
            //std::cout << "number of collisions: " << number_of_collisions << std::endl;
            if (pos_buffer != 0)
                cdslib::SaveValue(out, buffer, (pos_buffer + 63) / 64);
            delete[] buffer;
        }

        //! Returns true if the value stored at pos in the Hash contains the same
        // information than h_value
        bool
        verify(size_type original_h_pos, size_type h_pos, size_type h_value) {
            if (h_pos == 0)  //the root
                return false;
            size_type pos = h_pos; //position in the hash table
            size_type mod_P_value = 0, n_value = 0;
            if (H[h_pos] != h_value)
                return false;
            pos = (pos + M - D.getValue(h_pos)) % M; //update position in case depending of the collisions
            if (pos != original_h_pos)
                return false;
            return true; //pos == original_h_pos  and H[h_pos] == h_value
        }

        void
        grow_tables(std::ostream &out) {
            size_type read_size_bits, max_size_bytes, read_space;
            size_type cont = 0, len_p, h_pos;
            size_type node_value, current_node_pos = 0;
            if (pos_buffer != 0) { //empty the buffer to out
                cdslib::SaveValue(out, buffer, (pos_buffer + 63) / 64);
                pos_buffer = 0;
            }
            delete [] buffer;
            uint64_t end_pos = out.tellp(); //last position
            size_type M_prime = (size_type)std::min(MAX_M, 2 * M);  //in practice is better to move by 2 than by (1+factor)
            ghlz78<t_width, D_type> new_ds(MAX_M, M_prime, sigma, alpha, bits_d, factor, bits_L);
            std::ifstream f_in(out_file_name, std::ios::in | std::ios::binary);
            f_in.clear();
            f_in.seekg(0, f_in.beg); //make sure that we are at the start of the file
            h_pos = LoadValue<uint64_t>(f_in); //read the first variable
            uint64_t current_pos = f_in.tellg();
            value_type c;
            bool need_data = true;
            uint64_t *reading_buffer;
            uint64_t pos_reading_buffer = 0;
            max_size_bytes = bits_L * (cdslib::BUFFER_SIZE / (8 * bits_L));
            out.seekp(0, out.beg); //move to the start
            SaveValue(out, (int64_t) 0); //save space to store the start of the data structure
            while (cont <  number_of_phrases) {
                if (need_data) {
                    read_space = end_pos - current_pos;
                    if (read_space > max_size_bytes)
                        read_space = max_size_bytes;
                    read_size_bits = 8 * read_space;
                    reading_buffer = (uint64_t *) (cdslib::LoadValue<char>(f_in, read_space));
                    current_pos += read_space;
                    need_data = false;
                }
                h_pos = cdslib::get_field(reading_buffer, pos_reading_buffer, pos_reading_buffer + bits_L - 1);
                auto phrase = extract(h_pos);
                pos_reading_buffer += bits_L;
                len_p = phrase.size();
                current_node_pos = 0;
                while (len_p > 0) {
                    //now insert in the new tables
                    node_value = current_node_pos * sigma + (phrase[len_p - 1]);
                    current_node_pos = new_ds.insert(node_value, out);
                    -- len_p;
                }
                if (pos_reading_buffer == read_size_bits) {
                    delete[] reading_buffer;
                    pos_reading_buffer = 0;
                    need_data = true;
                }
                ++cont;
            }
            if (pos_reading_buffer > 0)
                delete[] reading_buffer;
            f_in.close();
            //swap data structures
            //size_type ram = size_in_bytes() + new_ds.size_in_bytes();
            this->swap(new_ds);
        }

        //!Extracts the value path from h_pos to the root into the out file
        void
        extract(size_type h_pos, std::ostream &out) {
            std::vector<value_type> phrase;
            phrase = extract(h_pos);
            value_type p;
            size_type i = phrase.size();
            max_length = (size_type)std::max(max_length, i);
            //write to the buffer now
            while (i > 0) {
                p = phrase[i - 1];
                if (p == 0) {
                   // std::cout << "h_pos: " <<  h_pos << "  ";
                   // std::cout << "Error phrase[i-1] == 0" << std::endl;
                }
                cdslib::set_field(buffer, pos_buffer, pos_buffer + t_width, alphabet[p - 1]);
                pos_buffer += t_width;
                cdslib::check_buffer(buffer, pos_buffer, out);
                -- i;
            }
        }

        size_type
        size_in_bytes() {
            size_type mem = 0;
            mem += 8 * sizeof(size_type); //parameters
            mem += sdsl::size_in_bytes(H);
            mem += sdsl::size_in_bytes(D);
            mem += sdsl::size_in_bytes(alphabet);
            //missing few things, but should not greatly change the final size for big files
            return mem;
        }

    };
}

#endif //CDSLIB_GHLZ78
