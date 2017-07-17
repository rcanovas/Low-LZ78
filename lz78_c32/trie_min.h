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

#define N_1 2
#define N_t 510 
#define ALPHA 0.95

#include "basics.h"
#include "heap.h"
#include "lztrie.h"


typedef struct strie_min {
   byte *trie;            // trie
   heap heaps[256];  
   uint nid;              // nr of nodes
   uint nbits;
} *trie_min;

typedef struct {
   int NPAR_P;
   int MSIZE_P;
   int NCH_P;
   int PAR_P;
   int FLAGS_P;
   int IDS_P;
   int LETTERS_P;
   int CHILDS_P;
   int Psize; // size of the page in bytes
   int maxsize; // maximum number of parenthese in the page 
} PagePos;


	// creates trie
trie_min createTrie(uint u, int N1, int Nt, float alpha, bool isrevtrie);
        // inserts word[0...] into pTrie and returns new text ptr
        // insertion proceeds until we get a new trie node
void insertTrie (trie_min pTrie, byte *word, int *m, FILE *fin, FILE *fout);
        // frees the trie
void destroyTrie (trie_min pTrie);
	// represents pTrie with parentheses, letters and ids
uint representTrie (trie_min pTrie, uint *parent, byte *letters, int *maxdepth, FILE *fp);

uint sizeTrie(trie_min pTrie, uint *maxdepth);
#endif
