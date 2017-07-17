/* Low-LZ78 - compressed data structures library
    Copyright (C) 2017 Diego Arroyuelo

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

// Implements operations over a sequence of balanced parentheses

#ifndef PARENTHESESINCLUDED
#define PARENTHESESINCLUDED

#include "basics.h"
#include "bitmap.h"
#include "hash.h"

typedef struct sparentheses
   { uint *data;  	// string
     bitmap bdata;   	// bitmap of string
     uint n;    	// # of parentheses
     uint sbits;    	// bits for near pointers
     hash sftable;	// table of far forward pointers
     hash sbtable;	// table of far backward pointers
     hash bftable;	// table of near forward pointers
     hash bbtable;	// table of near backward pointers
   } *parentheses;

	// creates a parentheses structure from a bitstring, which gets owned
        // n is the total number of parentheses, opening + closing
	// bwd says if you will want to perform findopen and enclose
parentheses createParentheses (uint *string, uint n, bool bwd);
	// frees parentheses structure, including the bitstream
void destroyParentheses (parentheses P);
	// the position of the closing parenthesis corresponding to (opening)
	// parenthesis at position i
uint findclose (parentheses P, uint i);
	// respectively, for closing parenthesis i
uint findopen (parentheses P, uint i);
	// # open - # close at position i, not included
uint excess (parentheses P, uint i);
	// open position of closest parentheses pair that contains the pair
	// that opens at i, ~0 if no parent
uint enclose (parentheses P, uint i);

uint sizeofParentheses(parentheses P);

#endif
