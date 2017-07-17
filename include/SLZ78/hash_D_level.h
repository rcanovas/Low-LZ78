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


// based on Cleary paper "Compact Hash Tables Using Bidirectional Linear Probing"

// Modified Poyias code (https://github.com/Poyias/mBonsai)
// Manage the D array using a hash table and a sublayer array

#ifndef CDSLIB_HASH_D_LEVEL_H
#define CDSLIB_HASH_D_LEVEL_H

#include "./../tools.h"
#include <limits.h>

namespace cdslib {

    class hash_D_level {
    public:
        typedef typename sdsl::int_vector<>::size_type  size_type;

    private:
        //sublayer displacement variables
        sdsl::int_vector<> D1; //hash table
        sdsl::int_vector<> V;
        sdsl::int_vector<> C;
        sdsl::int_vector<> satData;
        size_type emptyLoc;    // value that indicates in D1 that a cell has not been used
        size_type prime;        // prime number used for the hashtable D1
        size_type alpha_1;      // alpha value of D1;
        size_type Msl;          // size of D1
        size_type curEmptySlot; // auxiliary variable

        size_type load_limit;   //maximum number of elements
        size_type d_inserted;
        size_type M_size;
        hash_D_level *extra_level;


    public:

        // Empty constructor
        hash_D_level() {
            extra_level = nullptr;
        }

        //t_d_bits must be smaller than 7
        hash_D_level(size_type M, size_type _Msl = 0) {
            size_type cmax;
            double f = 0.062303;
            Msl = _Msl;
            if (Msl == 0)
                Msl = 1.2 * (((double)M * f) + 1.0);
            load_limit = (size_type) (std::ceil(Msl * (1.0 / (1.0 + f))));
            d_inserted = 0;
            M_size = M;
            extra_level = nullptr;
            cmax = M - 1;
            prime = cdslib::nearest_prime(cmax);
            alpha_1 = ULONG_MAX / prime; //taken from original implementation
            emptyLoc = (cmax / Msl) + 2;
            D1 = sdsl::int_vector<>(Msl, emptyLoc, sdsl::bits::hi(emptyLoc) + 1);
            satData = sdsl::int_vector<>(Msl, 0, 7); //using 7 bits as in the original implementation
            V = sdsl::int_vector<>(Msl, 0, 1);
            C = sdsl::int_vector<>(Msl, 1, 1);
        }


        void
        swap(hash_D_level& hdl) {
            if (this != &hdl) {
                D1.swap(hdl.D1); //hash table
                V.swap(hdl.V);
                C.swap(hdl.C);
                satData.swap(hdl.satData);
                std::swap(emptyLoc, hdl.emptyLoc);
                std::swap(prime, hdl.prime);
                std::swap(alpha_1, hdl.alpha_1);
                std::swap(Msl, hdl.Msl);
                std::swap(curEmptySlot, hdl.curEmptySlot);
                std::swap(load_limit, hdl.load_limit);
                std::swap(d_inserted, hdl.d_inserted);
                std::swap(M_size, hdl.M_size);
                std::swap(extra_level, hdl.extra_level);
            }
        }


        size_type
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0, e_bytes = 0;
            written_bytes += write_member(prime, out, child, "prime value for D1");
            written_bytes += write_member(alpha_1, out, child, "alpha value for D1");
            written_bytes += write_member(emptyLoc, out, child, "empty location value");
            written_bytes += write_member(Msl, out, child, "size D1");
            e_bytes += D1.serialize(out, child, "D1 hash table");
            e_bytes += V.serialize(out, child, "V bitarray");
            e_bytes += C.serialize(out, child, "C bitarray");
            e_bytes += satData.serialize(out, child, "satData array");
            bool need_more = false;
            if (extra_level != nullptr) {
                need_more = true;
                written_bytes += write_member(need_more, out, child, "extra table D1");
                written_bytes += (*extra_level).serialize(out, child, "Hash table");
            }
            else
                written_bytes += write_member(need_more, out, child, "extra table D1");
            written_bytes += e_bytes;
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void
        load(std::istream& in) {
            sdsl::read_member(prime, in);
            sdsl::read_member(alpha_1, in);
            sdsl::read_member(emptyLoc, in);
            sdsl::read_member(Msl, in);
            D1.load(in);
            V.load(in);
            C.load(in);
            satData.load(in);
            bool need_more = false;
            sdsl::read_member(need_more, in);
            if (need_more) {
                extra_level = new hash_D_level();
                extra_level->load(in);
            }
        }


        size_type
        getValue(size_type pos) {
            size_type res = find(pos);
            if (res < 135)
                return res;
            if(extra_level != nullptr)
                res = extra_level->getValue(pos);
            return res;
        }

        void
        setValue(size_type pos, size_type d) {
            if (d_inserted == load_limit) {
                if (extra_level == nullptr) {
                    //std::cout << "Needed another D1 level" << std::endl;
                    extra_level = new hash_D_level(M_size);
                }
                extra_level->setValue(pos, d);
            }
            else {
                insert_in_D1(pos, d);
                ++ d_inserted;
            }
        }


        bool
        is_full() {
            if (d_inserted == load_limit)
                return true;
            return false;
        }

        size_type
        get_size() {
            return Msl;
        }

    private:

        void
        insert_in_D1(size_type pos, size_type d) {
            size_type changeBit = 0;
            size_type cRand = cdslib::my_mod2(pos % prime, alpha_1, prime);
            size_type pos_hash = cRand % Msl;
            size_type value_hash = cRand / Msl;
            if (D1[pos_hash] == emptyLoc) {
                D1[pos_hash] = value_hash;
                satData[pos_hash] = d;
                V[pos_hash] = 1;
                C[pos_hash] = 1;
            }
            else {
                changeBit = getChangeBitLoc(pos_hash);
                if (V[pos_hash]==0) {  // not sure when this happen
                    startNewBlock(changeBit); // if changeBit == Msl this does nothing
                    V[pos_hash] = 1;
                    C[curEmptySlot] = 1;
                    D1[curEmptySlot] = value_hash;
                    satData[curEmptySlot] = d;
                }
                else {
                    findSpace(changeBit, value_hash);
                    D1[curEmptySlot] = value_hash;
                    C[curEmptySlot] = 0;
                    satData[curEmptySlot] = d;
                }
            }
        }

        size_type
        getChangeBitLoc(size_type curAddress) {
            size_type vOnesDown = 0, cOnesUp = 0;
            size_type count = 0;
            if (V[curAddress] == 1) // used position
                vOnesDown++;
            curAddress = (curAddress + Msl - 1) % Msl; //move on position downwards
            while(D1[curAddress] != emptyLoc) { // go downwards until empty slot
                count++;
                if (count >= D1.size()) {
                    std::cout << count << " D1 size: " << D1.size() << "  lf: " << load_limit << " d_i: "
                              << d_inserted << std::endl;
                    std::cout << "ERROR: table D1 is full"<< std::endl;
                    exit(0);
                }

                if (V[curAddress] == 1)
                    vOnesDown++;
                if (curAddress == 0) curAddress = Msl - 1;
                else                 curAddress--;
            }
            curEmptySlot = curAddress; // get empty slot
            if (vOnesDown == 0)
                return Msl; //value was not found
            curAddress = (curAddress + 1) % Msl; //move on position upward
            while (cOnesUp < vOnesDown) { // go upwards until conesUp==vOnes
                if (C[curAddress] == 1)
                    cOnesUp++;
                if (curAddress == Msl-1) curAddress = 0;
                else                     curAddress++;
            }
            //return associated C value
            return (curAddress + Msl - 1) % Msl;
        }

        // In case we want to start a new Block and the changeBit lies in another block.
        // We handle this situation differently.
        void
        startNewBlock(size_type changeBit) {
            if (changeBit == Msl)
                return;
            size_type tmpSlot, curC;
            curC = (changeBit + 1) % Msl;
            while (C[curC] == 0) {  //find the first position on C such that have value 1
                if (curC == Msl - 1) curC = 0;
                else                 curC++;
            }
            while (curEmptySlot != curC) { //swap values until arriving to curC
                if (curEmptySlot == Msl - 1) tmpSlot = 0;
                else                         tmpSlot = curEmptySlot + 1;
                D1[curEmptySlot] = D1[tmpSlot];
                C[curEmptySlot] = C[tmpSlot];
                satData[curEmptySlot] = satData[tmpSlot];
                curEmptySlot = tmpSlot;
            }
            curEmptySlot = (curEmptySlot + Msl - 1) % Msl;
        }

        bool
        itemExists(size_type cVal, size_type hash_value) {
            if (D1[cVal] == hash_value)
                return true;
            return false;
        }

        // We already found our associateC location.
        // We push c values and hashtables accordingly,
        // so we make room to insert the new node in the associateC location
        // return boolean if we found the item in the
        void
        findSpace(size_type cVal, size_type hash_value) {
            size_type curC, tmpSlot;
            if ((itemExists(cVal, hash_value)) || (D1[cVal] == emptyLoc)) { //check if the value is already inserted
                curEmptySlot = cVal;
                return;
            }
            // start going upwards until block ends where c!=0
            curC = (cVal + 1) % Msl;
            //go upwards towards the end of the block
            while (C[curC] == 0) {
                if (itemExists(cVal, hash_value)) {
                    curEmptySlot = curC;
                    return;
                }
                if (curC == Msl - 1)
                    curC = 0;
                else curC++;
            }
            if (curC == 0)
                curC = Msl - 1;
            else
                curC--;        //go one back to stay in the block
            while (curEmptySlot != curC) { // push all the slots up to curC to insert it in curC
                if (curEmptySlot == Msl - 1)
                    tmpSlot = 0;
                else
                    tmpSlot = curEmptySlot + 1;
                D1[curEmptySlot] = D1[tmpSlot];
                C[curEmptySlot] = C[tmpSlot];
                satData[curEmptySlot] = satData[tmpSlot];
                if (curEmptySlot == Msl - 1)
                    curEmptySlot = 0;
                else
                    curEmptySlot++;
            }
            curEmptySlot = curC;
        }

        size_type
        find(size_type pos) {
            size_type cRand = cdslib::my_mod2(pos % prime, alpha_1, prime);
            size_type pos_hash = cRand % Msl;
            size_type value_hash = cRand / Msl;

            if (V[pos_hash] == 0)
                return 135;  //not stored
            else {
                size_type exists = getSatelite(pos_hash, (getChangeBitLoc(pos_hash)), value_hash);
                if (exists != Msl) //found the value
                    return satData[curEmptySlot];
                return 135; //not stored here
            }
        }

        // Returns the location of the displacement value
        // which is stored as satelite data.
        // Called in find only when node exists.
        size_type
        getSatelite(size_type vVal, size_type cVal, size_type value_hash) {
            size_type curC, tmpSlot;
            if ((itemExists(cVal, value_hash))) { //check if the value is already inserted
                curEmptySlot = cVal;
                return cVal;
            }
            else if (D1[cVal] == emptyLoc) {
                curEmptySlot = cVal;
                return Msl;
            }
            curC = (cVal + 1) % Msl; // start going upwards until block ends where c!=0
            while (C[curC] == 0) { //go upwards towards the end of the block
                if (itemExists(curC, value_hash)) {
                    curEmptySlot = curC;
                    return curC;
                }
                if (curC == Msl - 1)
                    curC = 0;
                else
                    curC++;
            }
            return Msl;
        }

    };

}

#endif //CDSLIB_HASH_D_LEVEL_H


