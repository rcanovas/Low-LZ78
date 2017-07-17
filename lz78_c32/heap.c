
  // Heap of fixed size objects, fast to alloc and especially to totally free
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


#include "heap.h"

#define BLK 2048

#ifdef MEMCTRL
extern int account,currmem,newmem,maxmem;
#endif

        // Creates a heap of elements of size siz

heap createHeap (uint siz)

   { heap H;
     H = malloc (sizeof(struct sheap));
     H->siz = siz;
     H->cap = BLK/siz; if (H->cap == 0) H->cap = 1;
     H->blocks = NULL;
     H->free = NULL;
     H->first = H->cap;
#ifdef MEMCTRL
     H->totsize = 0;
#endif
     return H;
   }

        // Gets a new element from H

void *mallocHeap (heap H)

   { void *elem;
     void **b;
#ifdef MEMCTRL
     H->totsize++;
     newmem = currmem+H->siz;
     if (newmem > maxmem) maxmem = newmem;
     if (currmem/1024 != newmem/1024)
        printf ("Memory: %i Kb, maximum: %i Kb\n",newmem/1024,maxmem/1024);
     currmem = newmem;
#endif
     if (H->free)
	{ elem = (void*)H->free;
	  H->free = *H->free;
	  return elem;
	}
     if (H->first != H->cap)
	{ elem = (void*)(((uint)(H->blocks+1)) + H->first * H->siz);
	  H->first++;
	  return elem;
	}
#ifdef MEMCTRL
     account = 0;
#endif
     b = malloc (sizeof(void*) + H->cap * H->siz);
#ifdef MEMCTRL
     account = 1;
#endif
     b[0] = H->blocks;
     H->blocks = b;
     elem = (void*)(H->blocks+1);
     H->first = 1;
     return elem;
   }

        // Frees ptr from heap H

void freeHeap (heap H, void *ptr)

   { void **p;
     if (!ptr) return;
#ifdef MEMCTRL
     H->totsize--;
     newmem = currmem-H->siz;
     if (currmem/1024 != newmem/1024)
        printf ("Memory: %i Kb, maximum: %i Kb\n",newmem/1024,maxmem/1024);
     currmem = newmem;
#endif
     p = (void**)ptr;
     *p = H->free;
     H->free = p;
   }

        // Frees everything in heap H

void destroyHeap (heap H)

   { void **b;
#ifdef MEMCTRL
     account = 0;
#endif
     while ((b = H->blocks) != NULL)
	{ H->blocks = b[0];
	  free (b);
	}
#ifdef MEMCTRL
     account = 1;
     newmem = currmem-H->totsize*H->siz;
     if (currmem/1024 != newmem/1024)
        printf ("Memory: %i Kb, maximum: %i Kb\n",newmem/1024,maxmem/1024);
     currmem = newmem;
#endif
     free (H);
   }
