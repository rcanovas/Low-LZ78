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

#ifndef LZ78CINCLUDED
#define LZ78CINCLUDED

#define _FILE_OFFSET_BITS 64

#include "index_uc.c"



int main (int argc, char **argv)
 { 
    byte *text;
    uint i,level;
    unsigned long u;
    FILE *fin, *fout;
    struct stat sdata;
    lz78Parsing LZ78;
    char fnamext[1024];
    
    if (argc != 2) { 
       fprintf (stderr,"Usage: %s <file> \n indexes <file> and produces files "
                       "<file>.*,\n       you can then delete <file>\n",argv[0]);
       exit (1);
    }
    // load text
    printf ("Compressing %s...\n",argv[1]);
    if (stat(argv[1],&sdata) != 0) {
        fprintf (stderr,"Cannot stat file %s\n",argv[1]);
        exit(1);
    }
    u =  sdata.st_size;
    // allocates memory for the text buffer
    text = malloc(sizeof(byte)*BUFFER_TEXT_SIZE);
    fin = fopen(argv[1], "r");
    if (fin == NULL) { 
       fprintf (stderr,"Cannot open file %s\n",argv[1]);
       exit(1);
    }
    txtload(fin, text);  // reads part of the text in a buffer
    // build index
    LZ78 = buildLZTrie(text,fin,u);
    free(text); // frees the buffer of the text
    if (fclose(fin) != 0) { 
       fprintf (stderr,"Error: cannot close file %s\n",argv[1]);
       exit(1);
    }
    // save index
    sprintf (fnamext,"%s.lzt",argv[1]);
    fout = fopen (fnamext,"w");
    if (fout == NULL) {
       fprintf (stderr,"Error: cannot create file %s\n",fnamext);
       exit(1);
    }
     
    fwrite(&LZ78.z, sizeof(uint), 1, fout);
    fwrite(&LZ78.depth, sizeof(uint), 1, fout);
    fwrite(LZ78.letters, sizeof(byte), LZ78.z, fout);
    uint aux = ((unsigned long long)LZ78.z*bits(LZ78.z-1)+64-1)/64;
    fwrite(LZ78.parent, sizeof(unsigned long long), aux, fout);

    if (fclose(fout) != 0) {
       fprintf (stderr,"Error: cannot write file %s\n",fnamext);
       exit(1);
    }
    // free structures
    return 0;
 }
#endif
