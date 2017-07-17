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

// Basics

// #include "basics.h" included later to avoid macro recursion for malloc
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
	// Memory management

#ifdef MEMCTRL
int account = 1;
int currmem = 0;
int newmem;
int maxmem = 0;
#endif

void *Malloc (int n)

   { void *p;
     if (n == 0) return NULL;
#ifndef MEMCTRL
     p = (void*) malloc (n);
#else
     p = (void*) (malloc (n+sizeof(int))+sizeof(uint));
#endif
     if (p == NULL)
        { fprintf (stderr,"Could not allocate %i bytes\n",n);
          exit(1);
        }
#ifdef MEMCTRL
     *(((int*)p)-1) = n;
     if (account)
        { newmem = currmem+n;
          if (newmem > maxmem) maxmem = newmem;
          if (currmem/1024 != newmem/1024)
             printf ("Memory: %i Kb, maximum: %i Kb\n",newmem/1024,maxmem/1024);
          currmem = newmem;
	}
#endif
     return p;
   }

void Free (void *p)

   { 
#ifndef MEMCTRL
     if (p) free (p);
#else
     if (!p) return;
     if (account)
        { newmem = currmem - *(((int*)p)-1);
          free ((void*)(((int)p)-sizeof(int)));
          if (currmem/1024 != newmem/1024)
             printf ("Memory: %i Kb, maximum: %i Kb\n",newmem/1024,maxmem/1024);
          currmem = newmem;
	}
#endif
   }

void *Realloc (void *p, int n)

   { if (p == NULL) return Malloc (n);
     if (n == 0) { Free(p); return NULL; }
#ifndef MEMCTRL
     p = (void*) realloc (p,n);
#else
     if (account)
        newmem = currmem - *(((int*)p)-1);
     p = (void*) (realloc ((void*)(((int)p)-sizeof(int)),n+sizeof(int))+sizeof(int));
     *(((int*)p)-1) = n;
#endif
     if (p == NULL)
        { fprintf (stderr,"Could not allocate %i bytes\n",n);
          exit(1);
        }
#ifdef MEMCTRL
     if (account)
        { newmem = newmem + n;
          if (newmem > maxmem) maxmem = newmem;
          if (currmem/1024 != newmem/1024)
             printf ("Memory: %i Kb, maximum: %i Kb\n",newmem/1024,maxmem/1024);
          currmem = newmem;
	}
#endif
     return p;
   }

#include "basics.h"

        // bits needed to represent a number between 0 and n

uint bits (uint n)

   { uint b = 0;
     while (n)
	{ b++; n >>= 1; }
     return b;
   }

        // returns e[p..p+len-1], assuming len <= W
/*
uint bitget (uint *e, uint p, uint len)

   { uint answ;
     e += p >> bitsW; p &= (1<<bitsW)-1;
     answ = *e >> p;
     if (len == W)
          { if (p) answ |= (*(e+1)) << (W-p);
          }
     else { if (p+len > W) answ |= (*(e+1)) << (W-p);
            answ &= (1<<len)-1;
          }
     return answ;
   }
*/


uint bitget (uint *e, uint p, uint len)

   { register uint i=p>>5, j=p&0x1F;
     if (j+len <= W)
        return (e[i] << (W-j-len)) >> (W-len);
     else 
        return (e[i] >> j) | ((e[i+1] << (64-j-len)) >> (W-len)); 
   }


uint bitgetLong (unsigned long long *e, unsigned long long p, uint len)
   
   { unsigned long long i=p>>6, j=p&0x3F;
     if (j+len <= 64)
        return (e[i] << (64-(j+len))) >> (64-len);
     else
        return (uint)(e[i] >> j) | (uint)((e[i+1] << (64-(j+len-64))) >> (64-len));
   }    



  	// writes e[p..p+len-1] = s, len <= W




void bitput (register uint *e, register uint p,
             register uint len, register uint s)

  {             uint i = p >> bitsW, j = j=p-i*W;
                uint mask = ((j+len) < W ? ~0u << (j+len) : 0)
                        | ((W-j) < W ? ~0u >> (W-j) : 0);
                e[i] = (e[i] & mask) | s << j;
                if (j+len>W) {
                        mask = ((~0u) << (len+j-W));
                        e[i+1] = (e[i+1] & mask)| s >> (W-j);
                }


/*     e += p >> bitsW; p &= (1<<bitsW)-1;
     if (len == W)
          { *e |= (*e & ((1<<p)-1)) | (s << p);
            if (!p) return;
            e++;
            *e = (*e & ~((1<<p)-1)) | (s >> (W-p));
          }
     else { if (p+len <= W)
               { *e = (*e & ~(((1<<len)-1)<<p)) | (s << p);
                 return;
               }
            *e = (*e & ((1<<p)-1)) | (s << p);
            e++; len -= W-p;
            *e = (*e & ~((1<<len)-1)) | (s >> (W-p));
          }*/
   }


void bitputLong (unsigned long long *e, unsigned long long p,
             register uint len, register uint s)

   { 
     e += p >> 6; p &= (1<<6)-1;
     if (len == 64)
          { *e |= (*e & ((1<<p)-1)) | (s << p);
            if (!p) return;
            e++;
            *e = (*e & ~((1<<p)-1)) | (s >> 64-p);
          }
     else { if (p+len <= 64)
               { *e = (*e & ~((((unsigned long long)1<<len)-1)<<p)) | ((unsigned long long)s << p);
                 return;
               }
            *e = (unsigned long long)(*e & (((unsigned long long)1<<p)-1)) | ((unsigned long long)s << p);
            e++; len -= 64-p;
            *e = (unsigned long long)(*e & ~(((unsigned long long)1<<len)-1)) | ((unsigned long long)s >> (64-p));
          }
   }





int bitgetint(uint *e, uint p, uint len)
 {
    int answ;
    e += p >> bitsW; p &= (1<<bitsW)-1;
    answ = *e >> p;
    if (len == W) {
       if (p) answ |= (*(e+1)) << (W-p);
    }
    else {
       if (p+len > W) answ |= (*(e+1)) << (W-p);
       answ &= (1<<len)-1;
    }
    if ((answ & (1<<(len-1))) != 0) answ |= (~(int)0)<<len;
    return answ;
 }

// writes e[p..p+len-1] = s, len <= W

void bitputint(register int *e, register uint p,
               register uint len, register uint s)
 {
    e += p >> bitsW; p &= (1<<bitsW)-1;
    s = s & ~(~0<<len);
    if (len == W) {
       *e |= (*e & ((1<<p)-1)) | (s << p);
       if (!p) return;
       e++;
       *e = (*e & ~((1<<p)-1)) | (s >> (W-p));
    }
    else {
       if (p+len <= W) {
          *e = (*e & ~(((1<<len)-1)<<p)) | (s << p);
          return;
       }
       *e = (*e & ((1<<p)-1)) | (s << p);
       e++; len -= W-p;
       *e = (*e & ~((1<<len)-1)) | (s >> (W-p));
    }
 }

void txtload (FILE *fin, byte *text)
 { 
    bzero(text, BUFFER_TEXT_SIZE); 
    uint aux = fread(text, 1, BUFFER_TEXT_SIZE-1, fin);
 }


