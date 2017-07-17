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


#ifndef LZ78UINCLUDED
#define LZ78UINCLUDED

#define _FILE_OFFSET_BITS 64


#include "lztrie.h"
#include "index_min.c"

int main (int argc, char **argv)
 { 
    byte *text;
    uint i, j, k, d, dAux/*, *ids*/, nbits, n, node;
    uint u;
    size_t read;
    FILE *fin, *fout;
    lz78Parsing LZ78;
    char fnamext[1024];
    unsigned long long ids;
     
    if (argc != 2) { 
       fprintf (stderr,"Usage: %s <file> \n indexes <file> and produces files "
	    	        "<file>.*,\n       you can then delete <file>\n",argv[0]);
       exit(1);
    }
    // loads compressed text
    printf ("Decompressing %s...\n",argv[1]);
    sprintf (fnamext,"%s.lzt",argv[1]);
    fin = fopen (fnamext,"r");
    if (fin == NULL) { 
       fprintf (stderr,"Error: cannot cread from file %s\n",fnamext);
       exit(1);
    }
    read = fread(&LZ78.z, sizeof(uint), 1, fin);
    read = fread(&LZ78.depth, sizeof(uint), 1, fin);
    LZ78.letters = malloc(LZ78.z*sizeof(byte));
    read = fread(LZ78.letters, sizeof(byte), LZ78.z, fin);
    uint aux = ((unsigned long long)LZ78.z*bits(LZ78.z-1)+64-1)/64;
    LZ78.parent = malloc(aux*sizeof(unsigned long long));
    read = fread(LZ78.parent, sizeof(unsigned long long), aux, fin);
    
    if (fclose(fin) != 0) { 
       fprintf (stderr,"Error: cannot read file %s\n",fnamext);
       exit(1);
    }
    
    text = malloc(sizeof(byte)*(LZ78.depth+1));
    
    nbits = bits(LZ78.z-1);
    
    fout = fopen(argv[1], "w");
    if (fout == NULL) { 
       fprintf (stderr,"Error: cannot create file %s\n",argv[1]);
       exit(1);
    }
    text[LZ78.depth] = '\0';
    for (i=1; i < LZ78.z; i++) {
       k = LZ78.depth - 1;
       for (j=i; j > 0; j=bitgetLong(LZ78.parent, (unsigned long long)j*nbits,nbits)) {
          text[k--] = LZ78.letters[j];
       }
       fprintf(fout,"%s",&text[k+1]);
    }

    free((void*)text);

    if (fclose(fout) != 0) { 
       fprintf (stderr,"Error: cannot write file %s\n",fnamext);
       exit(1);
    }
    // free structures
    free(LZ78.letters);
    free(LZ78.parent);

    return 0;
 }
#endif
