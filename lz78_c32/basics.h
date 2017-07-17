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

#ifndef BASICSINCLUDED
#define BASICSINCLUDED


  // Includes 

#include <sys/types.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <unistd.h>

  // Memory management

#define malloc(n) Malloc(n)
#define free(p) Free(p)
#define realloc(p,n) Realloc(p,n)

void *Malloc (int n);
void Free (void *p);
void *Realloc (void *p, int n);

  // Data types

typedef unsigned char byte;
// typedef unsigned int uint;

typedef int bool;
#define true 1
#define false 0

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

  // Bitstream management

#define W (8*sizeof(uint))
#define bitsW 5 // OJO

	// bits needed to represent a number between 0 and n
uint bits (uint n);
        // returns e[p..p+len-1], assuming len <= W
uint bitget (uint *e, uint p, uint len);
        // writes e[p..p+len-1] = s, assuming len <= W
void bitput (uint *e, uint p, uint len, uint s);

uint bitgetLong (unsigned long long *e, unsigned long long p, uint len);
        // writes e[p..p+len-1] = s, assuming len <= W
void bitputLong (unsigned long long *e, unsigned long long p, uint len, uint s);

int bitgetint(uint *e, uint p, uint len);

void bitputint(register int *e, register uint p,
               register uint len, register uint s);
         // reads a part of the text in a buffer
void txtload (FILE *fin, byte *text);

         // reads bit p from e
#define bitget1(e,p) ((e)[(p)/W] & (1<<((p)%W)))
	// sets bit p in e
#define bitset(e,p) ((e)[(p)/W] |= (1<<((p)%W)))
	// cleans bit p in e
#define bitclean(e,p) ((e)[(p)/W] &= ~(1<<((p)%W)))

// size of the buffer for the text
#define BUFFER_TEXT_SIZE 1240000
  
#endif
