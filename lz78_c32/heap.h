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


  // Heap of fixed size objects, fast to alloc and especially to totally free

#ifndef HEAPINCLUDED
#define HEAPINCLUDED

#include "basics.h"

typedef struct sheap
   { uint siz;      // element size
     uint cap;      // block capacity in # of elems
     void **blocks; // list of blocks, next is block[0]
     void **free;   // list of free positions, each contains next
     uint first;    // first free position in current block
     uint totsize;  // total memory allocated here
   } *heap;

	// Creates a heap of elements of size siz
heap createHeap (uint siz);
	// Gets a new element from H
void *mallocHeap (heap H);
	// Frees ptr from heap H
void freeHeap (heap H, void *ptr);
	// Frees everything in heap H
void destroyHeap (heap H);

#endif
