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



#ifndef CDSLIB_TOOLS_H
#define CDSLIB_TOOLS_H


#include <sdsl/bit_vectors.hpp>
#include <map>
#include <iostream>
#include <fstream>


namespace cdslib {

    const uint64_t BUFFER_SIZE = 8388608; //1 megabytes in bits
                                          //83886080 -->10 megabytes in bits


    //! Prints the text length, its entropy and information of
    // the number of words in case that the case contains
    // a set of string separated by a given symbol
    void
    get_stats_text(sdsl::int_vector<> text, size_t separator = (size_t)-1) {
        double ent = 0;
        typedef typename sdsl::int_vector<>::size_type size_type;
        std::map<size_type, size_type> t_alphabet;
        size_type n = text.size(), n_alph = 0;
        size_type num_words = 0, length = 0;
        size_type occ = 0;
        for (size_type i = 0; i < n; ++i) {
            occ = t_alphabet[text[i]] + 1;
            t_alphabet[text[i]] = occ;
            if (text[i] == separator)
                ++num_words;
            else
                ++length;
        }
        n_alph = t_alphabet.size();
        for(auto it = t_alphabet.begin(); it != t_alphabet.end(); ++it)
            ent += ((it->second * 1.0) / (n * 1.0)) * log2((n * 1.0) / (it->second * 1.0));
        std::cout << "Text Length: " << text.size() << std::endl;
        std::cout << "H0: " << ent << std::endl;
        std::cout << "Number of words: " << num_words << "  Average size: " << (length * 1.0) / num_words  << std::endl;
    }

    //! Returns the nearest prime number p>n
    size_t
    nearest_prime(size_t n) {
        size_t prime = n + 1;
        size_t sqrt_prime, i;
        //we suppose that prime will not be longer than the maximum  unsigned int value
        while (true) {
            if (prime % 2 != 0) {
                sqrt_prime = (size_t) (std::sqrt(prime) + 1);
                for (i = 3; i < sqrt_prime; i += 2) {
                    if (prime % i == 0)
                        break;
                }
                if (i >= sqrt_prime)
                    break;
            }
            prime++;
        }
        return prime;
    }

    //! Returns the modular multiplicative inverse of alpha under modulo m
    size_t
    mod_mul_inv(size_t alpha, size_t m) {
        int64_t mod = m, x0 = 0, x1 = 1;
        int64_t tmp_swap, div;
        if (m == 1)
            return 1;
        while (alpha > 1) {
            div = alpha / m;
            tmp_swap = m;
            m = alpha % m;
            alpha = tmp_swap;
            tmp_swap = x0;
            x0 = x1 - div * x0;
            x1 = tmp_swap;
        }
        if (x1 < 0)
            x1 += mod;
        return x1;
    }

    //! Computes (x*y)%M in case that x*y > uint64_t.
    //We assume that 2 * M always is smaller than the maximum value
    //that can be stored into a uint64_t
    uint64_t
    my_mod(uint64_t x, uint64_t y, uint64_t M) {
        //std::cout << "1 ";
        if (x == 0 or y == 0)
            return 0;
        uint64_t z = (UINT64_MAX / x) - 1;
        if (z > y) {
            return (x * y) % M;
        }
        else {
            y = my_mod(x, y - z, M);
            x = (x * z) % M; //should be possible now
            return (x + y) % M;
        }
    }


     //! Computes (x*y)%M in case that x*y > uint64_t.
    //We assume that 2 * M always is smaller than the maximum value
    //that can be stored into a uint64_t
    uint64_t
    my_mod2(uint64_t x, uint64_t y, uint64_t M) {
         //std::cout << "1 ";
         if (x == 0 or y == 0)
             return 0;
         uint64_t z = (UINT64_MAX / x);
         if (z > y) {
             return (x * y) % M;
         }
         else {
             //--> xy == xz + x(y-z)
             uint64_t w = y / z;
             //--> xy== xz + x(wz + (y - wz) - z)--> xz+xzw+x(y-wz)-xz --> xzw + x(y-wz)
             y = (x * (y - z * w )) % M;
             x = (x * z) % M; //should be possible now
             x = my_mod2(x, w, M);
             return (x + y) % M;
         }
     }


    //! Store a given bitsequence into an uint64_t array A s
    // @param A Array
    // @param ini Starting position
    // @param fin End position (would write until fin)
    // @param x Value to be stored
    inline void
    set_field(uint64_t *A, const uint64_t ini, const uint64_t fin, const uint64_t x) {
        if (ini > fin)
            return;
        uint64_t i = ini / 64, j = ini - i * 64;
        uint64_t len = (fin - ini + 1);
        uint64_t mask = ((j + len) < 64 ? ~(uint64_t)0 << (j + len) : 0)
                      | ((64 - j) < 64 ? ~(uint64_t)0 >> (64 - j) : 0);
        A[i] = (A[i] & mask) | x << j;
        if (j + len > 64) {
            mask = ((~(uint64_t)0) << (len + j - 64));
            A[i + 1] = (A[i + 1] & mask) | x >> (64 - j);
        }
    }


    //! Retrieve a given index from array A
    // @param A Array
    // @param ini Starting position
    // @param fin End position (would read until fin -1)
    // @param index Position to be retrieved
    inline uint64_t
    get_field(const uint64_t *A, const uint64_t ini, const uint64_t fin) {
        if (ini > fin)
            return (uint64_t)-1;
        uint64_t i = ini / 64, j = ini - 64 * i;
        uint64_t result;
        uint64_t len = (fin - ini + 1);
        if (j + len <= 64)
            result = (A[i] << (64 - j - len)) >> (64 - len);
        else {
            result = A[i] >> j;
            result = result | (A[i + 1] << (128 - j - len)) >> (64 - len);
        }
        return result;
    }

    //! Saves len values into an ostream.
    template <typename T> void
    SaveValue(std::ostream &out, T *val, const uint64_t length) {
        assert(out.good());
        out.write(reinterpret_cast<char *>(val), length * sizeof(T));
    }

    //! Loads len values from an istream.
    template <typename T> T *
    LoadValue(std::istream &in, const uint64_t length) {
        assert(in.good());
        T *ret = new T[length];
        in.read(reinterpret_cast<char *>(ret), length * sizeof(T));
        return ret;
    }

    //! Saves a value into an ostream.
    template <typename T> void
    SaveValue(std::ostream &out, T val) {
        assert(out.good());
        out.write(reinterpret_cast<char *>(&val), sizeof(T));
    }

    //! Loads a value from an istream.
    template <typename T> T
    LoadValue(std::istream &in) {
        assert(in.good());
        T ret;
        in.read(reinterpret_cast<char *>(&ret), sizeof(T));
        return ret;
    }

    //! Initiallizes a buffer of aprox size BUFFER_SIZE
    inline uint64_t *
    init_buffer() {
        uint64_t  size = (uint64_t)( (1.2 * BUFFER_SIZE + 64) / 64);
        uint64_t *buffer = new uint64_t[size]; //ask for 20% in case of overflow
        for (uint64_t i = 0; i < size; ++i)
            buffer[i] = 0;
        return buffer;
    }


    inline void
    check_buffer(uint64_t *buffer, uint64_t  &pos_buffer, std::ostream &fp) {
        uint64_t copy = 0, pos = 0;
        if (pos_buffer >= BUFFER_SIZE) {
            copy = BUFFER_SIZE / 64;
            SaveValue(fp, buffer, copy);
            while ((copy + pos) * 64 < pos_buffer) {
                buffer[pos] = buffer[copy + pos];
                pos++;
            }
            pos_buffer = pos_buffer - copy * 64;
        }
    }

    uint64_t
    write_vector_selection(sdsl::int_vector<> &data, sdsl::int_vector<> &bitmap,
                           uint64_t n, std::ostream &out) {
        uint64_t *buffer = init_buffer();
        uint64_t pos_buffer = 0, bpv = data.width(), len = bitmap.size();
        SaveValue(out, n);
        SaveValue(out, bpv);
        for (uint64_t i = 0; i < len ; ++i) {
            if (bitmap[i]) {
                cdslib::set_field(buffer, pos_buffer, pos_buffer + bpv, data[i]);
                pos_buffer += bpv;
                cdslib::check_buffer(buffer, pos_buffer, out);
            }
        }
        if (pos_buffer != 0)
            SaveValue(out, buffer, (pos_buffer + 63) / 64);
        std::cout << "Space used by sel. values: " << 8 * ((n * bpv + 63) / 64) << " bytes" << std::endl;
        delete [] buffer;
        return 8 * ((n * bpv + 63) / 64);
    }

    void
    load_vector_selection(sdsl::int_vector<> &data, sdsl::int_vector<> &bitmap,
                          uint64_t default_value, std::istream &in) {
        uint64_t *buffer = init_buffer();
        uint64_t pos_buffer = 0, len = bitmap.size(), value = 0;
        uint64_t read_space = 0, read_size_bits = 0;
        uint64_t n = LoadValue<uint64_t>(in);
        uint64_t bpv = LoadValue<uint64_t>(in);
        sdsl::int_vector<> array(len, default_value, bpv);
        uint64_t max_size_bytes = bpv * (BUFFER_SIZE / (8 * bpv));
        uint64_t current_pos = in.tellg();
        uint64_t end_pos = current_pos + 8 * ((n * bpv + 63) / 64);
        bool need_data = true;
        for (uint64_t i = 0; i < len ; ++i) {
            if (bitmap[i]) {
                if (need_data) {
                    read_space = end_pos - current_pos;
                    if (read_space > max_size_bytes)
                        read_space = max_size_bytes;
                    read_size_bits = 8 * read_space;
                    buffer = (uint64_t *)(LoadValue<char>(in, read_space));
                    current_pos += read_space;
                    need_data = false;
                }
                value = get_field(buffer, pos_buffer, pos_buffer + bpv - 1);
                pos_buffer += bpv;
                if (pos_buffer == read_size_bits) {
                    delete [] buffer;
                    pos_buffer = 0;
                    need_data = true;
                }
                array[i] = value;
                --n;
                if (n == 0)
                    break;
            }
        }
        if (pos_buffer != 0)  //should never happen
            delete[] buffer;
        data.swap(array);
    }
}


#endif //CDSLIB_TOOLS_H
