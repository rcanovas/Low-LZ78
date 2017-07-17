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

// Implements the LZtrie data structure

#ifndef LZTRIEINCLUDED
#define LZTRIEINCLUDED

#include "basics.h"
#include "parentheses.h"

typedef uint trieNode; // a node of lztrie

#define NULLT ((trieNode)(~0)) // a null node
#define ROOT ((trieNode)(0)) // the root node

typedef struct slztrie
   { uint *data;	// bitmap data
     parentheses pdata; // parentheses structure
     uint n;    	// # of parentheses
     byte *letters;	// letters of the trie
     uint nbits;	// log n
     unsigned long long *id;		// ids of the trie
     trieNode *boost;   // direct pointers to first children;
     uint u;
     uint maxdepth;
   } *lztrie;

	// creates a lztrie structure from a parentheses bitstring,
	// a letter array in preorder, and an id array in preorder,
        // all of which except the latter become owned
        // n is the total number of nodes (n letters/ids, 2n parentheses)
lztrie createLZTrie (uint *string, byte *letters, unsigned long long *id, uint n, uint u, uint maxdepth);
	// frees LZTrie structure, including the owned data
void destroyLZTrie (lztrie T);
        // stores lztrie T on file f
void saveLZTrie (lztrie T, FILE *f);
        // loads lztrie T from file f
lztrie loadLZTrie (FILE *f);
	// letter by which node i descends
byte letterLZTrie (lztrie T, trieNode i);
	// go down by letter c, if possible
trieNode childLZTrie (lztrie T, trieNode i, byte c);
	// go up, if possible
trieNode parentLZTrie (lztrie T, trieNode i);
	// subtree size
uint subtreesizeLZTrie (lztrie T, trieNode i);
	// smallest rank in subtree
uint leftrankLZTrie (lztrie T, trieNode i);
	// largest rank in subtree
uint rightrankLZTrie (lztrie T, trieNode i);
	// id of node
uint idLZTrie (lztrie T, trieNode i);
	// rth of position
uint rthLZTrie (lztrie T, uint pos);
        // is node i ancestor of node j?
bool ancestorLZTrie (lztrie T, trieNode i, trieNode j);
	// next node from i, in preorder, adds/decrements depth accordingly
	// assumes it *can* go on!
trieNode nextLZTrie (lztrie T, trieNode i, uint *depth);

uint getTextLengthLZTrie(lztrie T);

uint maxdepth(lztrie T);

uint sizeofLZTrie(lztrie T);

#endif
