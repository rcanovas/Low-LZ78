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


 // LZ78 trie data structure

#ifndef TRIEINCLUDED
#define TRIEINCLUDED

#include "basics.h"
#include "heap.h"
#include <string.h>

typedef struct striebody_uc
   { uint id;	// node id
     short nchildren;
     struct schild
       { byte car;
         struct striebody_uc *trie;
       } *children;
   } triebody_uc;

typedef struct strie_uc
   { triebody_uc trie;	// trie
     heap heaps[257];	// heaps
     uint nid;    	// nr of nodes
   } *trie_uc;

	// creates trie
trie_uc createTrie (void);
        // inserts word[0...] into pTrie and returns new text ptr
        // insertion proceeds until we get a new trie node
void insertTrie (trie_uc pTrie, byte *word, uint *m, FILE *fin, FILE *fout);
        // frees the trie
void destroyTrie (trie_uc pTrie);

uint sizeTrie(triebody_uc *t, uint *maxdepth, uint mydepth);
#endif
