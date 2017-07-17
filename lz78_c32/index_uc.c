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

 // Indexing module
#ifndef INDEXINCLUDED
#define INDEXINCLUDED

#include "trie_uc.h"

typedef struct {
   uint z;
   uint depth;
   byte *letters;
   unsigned long long *parent;
} lz78Parsing;


	// creates lztrie over a null-terminated text
	// it also creates *ids

lz78Parsing buildLZTrie (byte *text, FILE *fin, long u) 
 { 
    trie_uc T;
    uint n, textpos = 0;
    uint *parent;
    byte *letters;
    uint sizeHLZT, maxdepth;
    FILE *fp;
    uint i, nbits, id, j, peak;
    byte *texto;
    unsigned long long aux, *ids;
    lz78Parsing LZ78;

    if ((fp = fopen("LZTIDSAUX__", "w"))==NULL) {
       printf("Error trying to write ids into disk...\n");
       exit(1);
    }

	// first creates a hrbp T for LZTrie
    T = createTrie();
    do { 
       if (textpos == (BUFFER_TEXT_SIZE-1)) {
          txtload(fin, text);
          textpos = 0;
       }    
       insertTrie(T,text,&textpos,fin,fp);
    }
    while (text[textpos-1]);
	// now we represent LZTrie from its hrbp
    n = T->nid; 
    
    maxdepth = 0;
    sizeHLZT = sizeTrie(&(T->trie), &maxdepth, 0);    

    printf(" ======================================================================================\n");
    
    printf("   > Number of LZ78 phrases: \t \t \t \t \t \t %u\n", n-1);
    
    printf("   > Step 1: Compact dynamic trie \t \t \t \t \t %u bytes\n", sizeHLZT);
    
    peak = sizeHLZT;

    destroyTrie(T);
    fclose(fp);
    fp = fopen("LZTIDSAUX__", "r");    

    letters = malloc(n*sizeof(byte));
    aux = ((unsigned long long)n*bits(n-1)+64-1)/64;
    ids = malloc((size_t)aux*sizeof(unsigned long long)); 

    ids[0] = 0;
    nbits = bits(n-1); 
    
    for (i = 1; i < n; i++) {
       aux = fread(&letters[i], sizeof(byte), 1, fp);
       aux = fread(&id, sizeof(uint), 1, fp);
       bitputLong(ids, (unsigned long long)i*nbits, nbits, id);
    }
     
    fclose(fp);
    remove("LZTIDSAUX__");
    aux = ((unsigned long long)n*bits(n-1)+64-1)/64;
    printf("   > Step 2: Letters + Phrase ids \t \t \t \t \t %u bytes\n", n*sizeof(byte) + (uint) aux*sizeof(unsigned long long));
    printf(" ======================================================================================\n");
    

    if (peak > n*sizeof(byte) + aux*sizeof(unsigned long long))
       printf("\n   >> The peak of memory usage is at Step 1\n");
    else
       printf("\n   >> The peak of memory usage is at Step 2\n");
    
    LZ78.z = n;
    LZ78.depth = maxdepth+1;
    LZ78.letters = letters;
    LZ78.parent = ids;
    return LZ78;
 }
#endif
