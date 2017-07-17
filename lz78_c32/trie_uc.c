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

#include "trie_uc.h"

	// creates trie

trie_uc createTrie (void)

   { trie_uc pTrie;
     uint i;
     pTrie = malloc (sizeof(struct strie_uc));
     pTrie->nid = 0;
     pTrie->trie.id = pTrie->nid++;
     pTrie->trie.nchildren = 0;
     pTrie->trie.children = NULL;
     pTrie->heaps[0] = createHeap(sizeof(triebody_uc));
     for (i=1;i<=256;i++)
	pTrie->heaps[i] = createHeap(i*sizeof(struct schild));
     return pTrie;
   }

	// frees the trie

void destroyTrie (trie_uc pTrie)

   { uint i;
     for (i=0;i<=256;i++) destroyHeap (pTrie->heaps[i]);
     free (pTrie);
   }

	// inserts word[0...] into pTrie and returns new text ptr
	// insertion proceeds until we get a new trie node

void insertTrie (trie_uc pTrie, byte *word, uint *m, FILE *fin, FILE *fout)
 { 
    triebody_uc *t = &pTrie->trie;
    triebody_uc *nt;
    struct schild *newc;
    int i,j;
    // traverse pTrie with word[0...]
    while (true) { 
       i = 0;
       while (i < t->nchildren) {
          if (t->children[i].car >= word[*m]) break;
          i++;
       }
       if ((i == t->nchildren) || (t->children[i].car > word[*m]))
          break;  // not found, get out
       t = t->children[i].trie;
       (*m)++;
       if ((*m) == (BUFFER_TEXT_SIZE-1)) {
          txtload(fin, word); 
	  (*m) = 0;
       }
    }
    // at this point we fell off the trie, which is guaranteed to occur
    // since the text finishes with the unique character 0
    newc = mallocHeap(pTrie->heaps[t->nchildren+1]);
    memcpy (newc,t->children,i*sizeof(struct schild));
    memcpy (newc+i+1,t->children+i,(t->nchildren-i)*sizeof(struct schild));
    freeHeap (pTrie->heaps[t->nchildren],t->children);
    
    fwrite (&word[*m], sizeof(byte), 1, fout);
    fwrite (&t->id, sizeof(uint), 1, fout);
    
    t->children = newc;
    t->children[i].car = word[*m];
    (*m)++;
    nt = mallocHeap (pTrie->heaps[0]);
    t->children[i].trie = nt;
    t->nchildren++;
    // new node created	
    nt->id = pTrie->nid++;
    nt->nchildren = 0;
    nt->children = NULL;
 }


uint sizeTrie(triebody_uc *t, uint *maxdepth, uint mydepth)
 {
    uint i, nbytes;

    nbytes = sizeof(triebody_uc) + t->nchildren*(sizeof(byte)+sizeof(triebody_uc *));
    
    // traverse children
    for (i=0;i<t->nchildren;i++) {
       nbytes += sizeTrie(t->children[i].trie, maxdepth, mydepth+1);
       if (mydepth > (*maxdepth)) *maxdepth = mydepth;
    }
    return nbytes;
 }


