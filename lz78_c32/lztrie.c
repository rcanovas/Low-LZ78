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

#include "lztrie.h"

	// creates a lztrie structure from a parentheses bitstring,
	// a letter array in preorder, and an id array in preorder,
	// all of which except the latter become owned
        // n is the total number of nodes (n letters/ids, 2n parentheses)

lztrie createLZTrie (uint *string, byte *letters, unsigned long long *id, uint n, uint u, uint maxdepth)

   { lztrie T;
     uint i;
     T = malloc (sizeof(struct slztrie));
     T->data = string;
     T->pdata = createParentheses (string,2*n,true);
     T->n = n;
     T->nbits = bits(n-1);
     T->letters = letters;
     T->id = id;
     T-> u = u; // text length
     T->maxdepth = maxdepth; // longest LZ78 phrase
	// boost first access
     T->boost = malloc (256*sizeof(trieNode));
     for (i=0;i<256;i++) T->boost[i] = NULLT;
     i = 1; // shortcut for first child of root
     while (i != 2*n-1) // shortcut for its closing parenthesis
        { T->boost[T->letters[i-rank(T->pdata->bdata,i)]] = i;
                // shortcut for leftrankLZTrie
          i = findclose(T->pdata,i)+1;
        }
     return T;
   }

	// frees LZTrie structure, including the owned data

void destroyLZTrie (lztrie T)

   { 
     destroyParentheses (T->pdata);
     free (T->id);
     free (T->letters);
     free (T);
   }

	// stores lztrie T on file f

void saveLZTrie (lztrie T, FILE *f)

   { unsigned long long aux;

     if (fwrite (&T->n,sizeof(uint),1,f) != 1)
	{ fprintf (stderr,"Error: Cannot write LZTrie on file\n");
	  exit(1);
	}
     if (fwrite (T->data,sizeof(uint),(2*T->n+W-1)/W,f) != (2*T->n+W-1)/W)
	{ fprintf (stderr,"Error: Cannot write LZTrie on file\n");
	  exit(1);
	}
     if (fwrite (T->letters,sizeof(byte),T->n,f) != T->n)
	{ fprintf (stderr,"Error: Cannot write LZTrie on file\n");
	  exit(1);
	}
     aux = ((unsigned long long)T->n*T->nbits+64-1)/64;
     if (fwrite (T->id,sizeof(unsigned long long),(size_t)aux,f) != aux)
	{ fprintf (stderr,"Error: Cannot write LZTrie on file\n");
	  exit(1);
	}
     if (fwrite (&T->u,sizeof(uint),1,f) != 1)
	{ fprintf (stderr,"Error: Cannot write LZTrie on file\n");
	  exit(1);
	}
     if (fwrite (&T->maxdepth,sizeof(uint),1,f) != 1)
        { fprintf (stderr,"Error: Cannot write LZTrie on file\n");
          exit(1);
        }
   }

	// loads lztrie T from file f

lztrie loadLZTrie (FILE *f)

   { lztrie T;
     uint i;
     unsigned long long aux;
     T = malloc (sizeof(struct slztrie));
     if (fread (&T->n,sizeof(uint),1,f) != 1)
	{ fprintf (stderr,"Error: Cannot read LZTrie from file\n");
	  exit(1);
	}
     T->data = malloc (((2*T->n+W-1)/W)*sizeof(uint));
     if (fread (T->data,sizeof(uint),(2*T->n+W-1)/W,f) != (2*T->n+W-1)/W)
	{ fprintf (stderr,"Error: Cannot read LZTrie from file\n");
	  exit(1);
	}
     T->pdata = createParentheses (T->data,2*T->n,true);
     T->nbits = bits(T->n-1);
     T->letters = malloc (T->n*sizeof(byte));
     if (fread (T->letters,sizeof(byte),T->n,f) != T->n)
	{ fprintf (stderr,"Error: Cannot read LZTrie from file\n");
	  exit(1);
	}
     aux = ((unsigned long long)T->n*T->nbits+64-1)/64;
     T->id = malloc ((size_t)aux*sizeof(unsigned long long));
     if (fread (T->id,sizeof(unsigned long long),(size_t)aux,f) != aux)
	{ fprintf (stderr,"Error: Cannot read LZTrie from file\n");
	  exit(1);
	}
     if (fread (&T->u,sizeof(uint),1,f) != 1)
	{ fprintf (stderr,"Error: Cannot read LZTrie from file\n");
	  exit(1);
	}
     if (fread (&T->maxdepth,sizeof(uint),1,f) != 1)
        { fprintf (stderr,"Error: Cannot read LZTrie from file\n");
          exit(1);
        }
	// boost first access
     T->boost = malloc (256*sizeof(trieNode));
     for (i=0;i<256;i++) T->boost[i] = NULLT;
     i = 1; // shortcut for first child of root
     while (i != 2*T->n-1) // shortcut for its closing parenthesis
        { T->boost[T->letters[i-rank(T->pdata->bdata,i)]] = i;
                // shortcut for leftrankLZTrie
          i = findclose(T->pdata,i)+1;
        }
     return T;
   }

uint getTextLengthLZTrie(lztrie T) 
 {
    return T->u;
 }

uint maxdepth(lztrie T) 
 {
    return T->maxdepth;
 }


        // letter by which node i descends
byte letterLZTrie (lztrie T, trieNode i)
   { return T->letters[i-rank(T->pdata->bdata,i)]; // shortcut leftrankLZTrie
   }
	
       // go down by letter c, if possible
#ifdef QUERYREPORT
int JUMPS = 0; // accesses to nodes to run childLZTrie
int CALLS = 0; // calls to childLZTrie
#endif

trieNode childLZTrie (lztrie T, trieNode i, byte c)

   { trieNode j;
     byte tc;
#ifdef QUERYREPORT
     CALLS++;
#endif
     if (i == ROOT) return T->boost[c];
     j = i+1;
     while (!bitget1(T->data,j)) // there is a child here
	{ tc = T->letters[j-rank(T->pdata->bdata,j)];
		// shortcut for leftrankLZTrie: is a child by letter c?
	  if (tc > c) break;
	  if (tc == c) return j;
	  j = findclose(T->pdata,j)+1;
#ifdef QUERYREPORT
          JUMPS++;
#endif
	}
     return NULLT; // no child to go down by c
   }

	// go up, if possible

trieNode parentLZTrie (lztrie T, trieNode i)

   { if (i == ROOT) return NULLT; // parent of root
     if (T->boost[T->letters[i-rank(T->pdata->bdata,i)]] == i) return ROOT;
     return enclose (T->pdata,i);
   }

	// subtree size

uint subtreesizeLZTrie (lztrie T, trieNode i)

   { trieNode j;
     j = findclose (T->pdata,i);
     return (j-i+1)/2;
   }

	// depth

uint depthLZTrie (lztrie T, trieNode i)

   { return excess (T->pdata,i);
   }
   
	// smallest rank in subtree

uint leftrankLZTrie (lztrie T, trieNode i)

   { return i-rank(T->pdata->bdata,i);
   }

	// largest rank in subtree

uint rightrankLZTrie (lztrie T, trieNode i)

   { trieNode j;
     j = findclose (T->pdata,i);
     return j-rank(T->pdata->bdata,j)-1;
   }

	// id of node

uint idLZTrie (lztrie T, trieNode i)

   { uint pos = leftrankLZTrie (T,i);
     return rthLZTrie (T,pos);
   }

	// rth of position

uint rthLZTrie (lztrie T, uint pos)
   { return bitget ((uint *)(T->id),(unsigned long long)pos*T->nbits,T->nbits);
   }

        // is node i ancestor of node j?

bool ancestorLZTrie (lztrie T, trieNode i, trieNode j)

   { return (i <= j) && (j < findclose (T->pdata,i));
   }

        // next node from i, in preorder, adds/decrements depth accordingly
	// assumes it *can* go on!

trieNode nextLZTrie (lztrie T, trieNode i, uint *depth)

   { uint *data = T->data;
     i++;
     while (bitget1(data,i)) { i++; (*depth)--; }
     (*depth)++;
     return i;
   }

uint sizeofLZTrie(lztrie T)
 {
    return sizeof(struct slztrie) +
           sizeofParentheses(T->pdata) +
           T->n*sizeof(byte) + // letters
           ((uint)((unsigned long long)T->n*T->nbits+64-1)/64)*sizeof(unsigned long long); // ids
 }

