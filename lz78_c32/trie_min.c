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

#include "trie_min.h"
#include <math.h>
#include <string.h>


byte T[256];  
PagePos N[1000]; 
int lastN;   // indica la pos en N correspondiente a las paginas
             // de tamanio maximo posible (N_t)
int F[512];

void precomputeT(void) 
 {
    int i, j, cant; register byte mb;
    for (i = 0; i < 256; i++) {
       mb = 1; cant = 0;
       for (j = 0; j < 8; j++) {
          cant+= ((mb & i) != 0);
	  mb <<= 1;
       }
       T[i] = cant;
    }
 }


// creates trie
trie_min createTrie(uint u, int N1, int Nt, float alpha, bool isrevtrie)
 {
 // npages: numero distinto de tamanios de paginas
 // N1: size (number of parentheses) of the smallest pages
 // nt: size (number of parentheses) of the greatest pages
 
    trie_min pTrie;
    uint i, j;
    precomputeT();  
    
    pTrie = malloc(sizeof(struct strie_min));
    pTrie->nid = 1;
    pTrie->nbits = isrevtrie+bits(u-1);
    
    N[0].maxsize  = N1;
    N[0].NPAR_P    = 0;
    N[0].MSIZE_P   = 1;
    N[0].NCH_P     = 2;
    N[0].PAR_P     = 3;
    N[0].FLAGS_P   = N[0].PAR_P + ceil((float)N[0].maxsize/8);
    N[0].IDS_P     = N[0].FLAGS_P + ceil(ceil((float)N[0].maxsize/8)/2);
    N[0].LETTERS_P = N[0].IDS_P + (ceil((pTrie->nbits * (float)(N[0].maxsize/2)/32))*32)/8;
    N[0].CHILDS_P  = N[0].LETTERS_P + N[0].maxsize/2;
    N[0].Psize     = N[0].CHILDS_P + sizeof(byte **);
    for(i = 1; N[i-1].maxsize < Nt; i++) {
       N[i].maxsize   = ceil((N[i-1].maxsize/alpha)/2)*2;  // para que siempre sea un numero par
       N[i].NPAR_P    = 0;
       N[i].MSIZE_P   = 1;
       N[i].NCH_P     = 2;
       N[i].PAR_P     = 3;
       N[i].FLAGS_P   = N[i].PAR_P + ceil((float)N[i].maxsize/8);
       N[i].IDS_P     = N[i].FLAGS_P + ceil(ceil((float)N[i].maxsize/8)/2);
       N[i].LETTERS_P = N[i].IDS_P + (ceil((pTrie->nbits * (float)(N[i].maxsize/2)/32))*32)/8;
       N[i].CHILDS_P  = N[i].LETTERS_P + N[i].maxsize/2;
       N[i].Psize     = N[i].CHILDS_P + sizeof(byte **);
    }
    lastN = i - 1;
    N[lastN].maxsize = Nt;
    N[lastN].NPAR_P    = 0;
    N[lastN].MSIZE_P   = 1;
    N[lastN].NCH_P     = 2;
    N[lastN].PAR_P     = 3;
    N[lastN].FLAGS_P   = N[lastN].PAR_P + ceil((float)N[lastN].maxsize/8);
    N[lastN].IDS_P     = N[lastN].FLAGS_P + ceil(ceil((float)N[lastN].maxsize/8)/2);
    N[lastN].LETTERS_P = N[lastN].IDS_P + (ceil((pTrie->nbits * (float)(N[lastN].maxsize/2)/32))*32)/8;
    N[lastN].CHILDS_P  = N[lastN].LETTERS_P + N[lastN].maxsize/2;
    N[lastN].Psize     = N[lastN].CHILDS_P + sizeof(byte **); 
    // computes the Fit mapping
    for (i = 0; i < Nt; i++) {
       for (j = 0; j < lastN;j++)
          if (N[j].maxsize >= i) break;
       F[i] = j;
    }    
    // creates the root page of size N1
    pTrie->trie = calloc(N[0].Psize, sizeof(byte));
    pTrie->trie[N[0].MSIZE_P] = 0;
    pTrie->heaps[0] = createHeap(1);
    for (i = 1; i < 256; i++)
       pTrie->heaps[i] = createHeap(i*sizeof(byte *));
    return pTrie;
 }

// frees the trie
void destroyTrie (trie_min pTrie) 
 {
    uint i;
    for (i = 0; i < 256; i++) destroyHeap(pTrie->heaps[i]);	      
    free(pTrie);
 }

#define setvars(i)\
 {\
    NPAR_P = N[i].NPAR_P;\
    MSIZE_P        = N[i].MSIZE_P;\
    NCH_P      = N[i].NCH_P;\
    PAR_P  = N[i].PAR_P;\
    FL_P        = N[i].FLAGS_P;\
    I_P          = N[i].IDS_P;\
    LET_P      = N[i].LETTERS_P;\
    CH_P       = N[i].CHILDS_P;\
    PAGE_SIZE      = N[i].Psize;\
 }

int select_subtrie(byte *t, byte *tnew, uint flaginsert, 
                   byte **tinsert, uint *posinsert, uint nbits) 
/*
t: bloque donde ocurre el overflow
tnew: nuevo bloque creado para resolver el overflow
flaginsert: flag de t correspondiente a la posicion de insercion del nuevo nodo 
tinsert: parametro de salida, apuntara al bloque donde debe realizarse la insercion luego de resolver el overflow
posinsert: bit dentro de tinsert donde debe ocurrir la insercion del nuevo nodo
*/
 {
    register byte b, mb, mbnew; 
    uint curbyte, curbit, letterpos, curflag, bitinit, parentflag;
    uint curbytenew, letterposnew, li, subtriesize, cant, i, npar,
         ids, idsnew,idinit;
    int excess;
    uint NPAR_P, MSIZE_P, NCH_P, PAR_P, FL_P, I_P, LET_P, CH_P, PAGE_SIZE; 
    // INVARIANTE: t y tnew son del mismo tamanio maximo, es decir
    // t[1] == tnew[1]
    setvars(t[1] /*da lo mismo tnew[1]*/);
    curbyte = PAR_P; mb = 2; curbit = 1; 
    npar = ((uint)t[NPAR_P])<<1;
    letterpos = curflag = 0;
    /* selecciono el subtrie a bajar */
    while (curbit < npar) {
       if ((mb & t[curbyte]) == 0) break;
       letterpos++; // next letter
       mb <<= 2; curbit += 2; curflag++;
       if (mb == 0) {mb = 2; curbyte++;}
    }
    
    if (curbit >= npar) { // quiere decir que no hay subtrie para bajar, 
                          // tengo que bajar el nuevo nodo
       *posinsert = 0; *tinsert = tnew;
       curbyte = (flaginsert>>3);
       t[FL_P + curbyte] |= (1 << (flaginsert%8));
       cant = 0;
       for (i = 0; i < curbyte; i++) cant += T[t[FL_P+i]];
       mb = 1;
       for (i = 0; i < (parentflag%8); i++)
          cant += ((mb & t[FL_P + curbyte]) != 0); 
       return cant;
    }
    /* ya encontro el subtrie, estoy parado en su primer '(' */
    b = 0; excess = 0; bitinit = curbit; parentflag = curflag++;
    curbytenew = PAR_P; letterposnew = LET_P;
    idsnew = 0; idinit= ids = letterpos+1;
    li = letterpos = LET_P + letterpos + 1; // 1st letter to copy to tnew
    subtriesize = 0; mbnew = 1; 
    while (excess != -1) {
       if ((mb & t[curbyte]) == 0) {
          tnew[letterposnew] = t[letterpos];
	  bitput(((uint*)&tnew[I_P]),idsnew*nbits,nbits,
	         bitget(((uint*)&t[I_P]),ids*nbits,nbits));
	  letterposnew++;
	  letterpos++; ids++; idsnew++;
	  excess++;
       }
       else {
          if (excess) {
	     b |= mbnew;
	     t[curbyte] &= ~mb;
	  }
	  excess--;
       }
       mb <<= 1; curbit++; subtriesize++; mbnew <<= 1;
       if (mb == 0) {mb = 1; curbyte++;}
       if (mbnew == 0) {
          mbnew = 1;
	  tnew[curbytenew] = b;
	  curbytenew++; b = 0;
       }
    }
    subtriesize--; curbit--;
    if (b != 0) tnew[curbytenew] = b;
    tnew[NPAR_P] = (subtriesize >> 1);
    /* borro en t el subtrie copiado a tnew */
    mb = 1 << (bitinit&0x07); curbyte = PAR_P + (bitinit >> 3);
    mbnew = 1 << (curbit&0x07); curbytenew = PAR_P + (curbit >> 3);
    for (i = curbit; i < npar; i++) {
       if ((mbnew & t[curbytenew]) == 0) {
          t[curbyte] &= ~mb;
	  t[li++] = t[letterpos++];
          bitput(((uint*)&t[I_P]),idinit*nbits,nbits,
                 bitget(((uint*)&t[I_P]),ids*nbits,nbits));
          ids++; idinit++;
       }
       else {
          t[curbyte] |= mb;
	  t[curbytenew] &= ~mbnew;
       }
       mb <<= 1; mbnew <<= 1;
       if (mb == 0) {mb = 1; curbyte++;}
       if (mbnew == 0) {mbnew = 1; curbytenew++;}
    } 
    /* elimino en t los flags corresp. al subtrie que copie a tnew */
    mb = 1 << (curflag&0x07); curbyte = FL_P + (curflag >> 3);
    mbnew = 1 << ((curflag+(subtriesize>>1))&0x07); 
    curbytenew = FL_P + ((curflag+(subtriesize>>1))>> 3);
    npar = t[NPAR_P];
    for (i = curflag+(subtriesize>>1); i < npar; i++) {
       if ((mbnew & t[curbytenew]) == 0)
          t[curbyte] &= ~mb;
       else {
          t[curbyte] |= mb;
	  t[curbytenew] &= ~mbnew;
       }
       mb <<= 1; mbnew <<= 1;
       if (mb == 0) {mb = 1; curbyte++;}
       if (mbnew == 0) {mbnew = 1; curbytenew++;}      
    }
    t[NPAR_P] -= (subtriesize>>1); // updates the number of nodes in t
    /* agrego un nuevo flag en 1 */
    curbyte = FL_P + (parentflag>>3);
    t[curbyte] |= (1 <<(parentflag & 0x07));
    /* calculo la posicion en donde insertar el puntero al nuevo bloque */
    cant = 0;
    for (i = FL_P; i < curbyte; i++) cant += T[t[i]];
    cant += T[(t[curbyte]&~(~0<<(parentflag&0x07)))];
    /* determino la nueva posicion en que debe insertarse 
       el nuevo nodo, y el bloque en donde hacerlo  */
    if ((bitinit + subtriesize) < (*posinsert)) {
       *posinsert = *posinsert - subtriesize;
       *tinsert = t; // insercion ocurre en subarbol original, en pos 
                     // posterior a la del subtrie que baje
    }
    else if (*posinsert < bitinit) 
       *tinsert = t; // pos de insercion es misma original, en bloque original
    else {
       *posinsert = *posinsert - bitinit;
       *tinsert = tnew; // insercion ocurre en subtrie que baje a nuevo bloque
    }
    /* retorno la posicion donde debe insertarse el nuevo bloque  */ 
    return cant;
 }



byte *grow(byte *t, byte **tparent)
 {
    byte *tnew;
    register uint i, j, aux, auxt, auxtnew;
    tnew = calloc(N[t[1]+1].Psize, sizeof(byte));
    tnew[1]                 = t[1] + 1;
    auxt                    = t[1];
    auxtnew                 = tnew[1];
    tnew[N[auxtnew].NPAR_P] = t[N[auxt].NPAR_P];
    tnew[N[auxtnew].NCH_P]  = t[N[auxt].NCH_P];
    /* copia los parenthesis */
    j   = N[auxtnew].PAR_P;
    i   = N[auxt].PAR_P;
    aux = N[auxt].FLAGS_P;
    for (; i < aux; i++) tnew[j++] = t[i];
    /* copia los flags */
    j   = N[auxtnew].FLAGS_P;
    i   = N[auxt].FLAGS_P; 
    aux = N[auxt].IDS_P;
    for (; i < aux; i++) tnew[j++] = t[i];
    /* copia los ids */
    j   = N[auxtnew].IDS_P;
    i   = N[auxt].IDS_P;
    aux = N[auxt].LETTERS_P;
    for (; i < aux; i++) tnew[j++] = t[i];
    /* copia las letras */
    j   = N[auxtnew].LETTERS_P;
    i   = N[auxt].LETTERS_P;
    aux = N[auxt].CHILDS_P;
    for (; i < aux; i++) tnew[j++] = t[i];
    /* copia los punteros a los hijos */
    tnew[N[auxtnew].CHILDS_P]   = t[N[auxt].CHILDS_P];
    tnew[N[auxtnew].CHILDS_P+1] = t[N[auxt].CHILDS_P+1];
    tnew[N[auxtnew].CHILDS_P+2] = t[N[auxt].CHILDS_P+2];
    tnew[N[auxtnew].CHILDS_P+3] = t[N[auxt].CHILDS_P+3];
    // cambia en el padre el puntero a t por un puntero a tnew
    *tparent = tnew;
    /* libera la memoria de la pagina vieja */
    free((void *)t);
    return tnew; /* retorna un puntero a la nueva pagina*/
 }

byte *shrink(byte *t, byte **tinsert)
 {
    byte *tnew;
    uint npar, k;
    register uint i, j, auxtnew, auxt;
    npar = ((uint)t[0])<<1;
    if (F[npar]!=t[1]) { // solo ajusta si es necesario
       /* pide nueva pagina */
       auxtnew = F[npar];
       tnew = calloc(N[auxtnew].Psize, sizeof(byte));
       /* asigna tamanio correcto */
       tnew[1] = auxtnew;
       auxt    = t[1];  
       /* copia numero de parentesis en la pagina */
       tnew[N[auxtnew].NPAR_P] = t[N[auxt].NPAR_P];
       /* copia numero de hijos en la pagina */
       tnew[N[auxtnew].NCH_P] = t[N[auxt].NCH_P];
       /* copia parenthesis */
       i = N[auxt].PAR_P;
       j = N[auxtnew].PAR_P;
       for (; j < N[auxtnew].FLAGS_P; j++) tnew[j] = t[i++];
       /* copia los flags */
       i = N[auxt].FLAGS_P;
       j = N[auxtnew].FLAGS_P;
       for (; j < N[auxtnew].IDS_P; j++) tnew[j] = t[i++];
       /* copia los ids */
       i = N[auxt].IDS_P;
       j = N[auxtnew].IDS_P;
       for (; j < N[auxtnew].LETTERS_P; j++) tnew[j] = t[i++];
       /* copia las letras */
       i = N[auxt].LETTERS_P;
       j = N[auxtnew].LETTERS_P;
       for (; j < N[auxtnew].CHILDS_P; j++) tnew[j] = t[i++];
       /* copia los punteros a los hijos */
       tnew[N[auxtnew].CHILDS_P]   = t[N[auxt].CHILDS_P];
       tnew[N[auxtnew].CHILDS_P+1] = t[N[auxt].CHILDS_P+1];
       tnew[N[auxtnew].CHILDS_P+2] = t[N[auxt].CHILDS_P+2];
       tnew[N[auxtnew].CHILDS_P+3] = t[N[auxt].CHILDS_P+3];
    }
    else tnew = t;
    if (*tinsert == t) *tinsert = tnew;
    if (tnew != t) free((void *)t); // libera la vieja pagina
    return tnew; // retorna puntero a la nueva pagina
 }


void insert_page(trie_min pTrie, byte *t, byte **tparent, uint bytepos, 
                 uint curbit, uint letterpos, uint curflag, int id, 
		 byte letter, bool isrevtrie)
 {
    uint npar, lastflagbyte, curflagbyte, cantflags, i, maxsize, 
         nbits, j;
    uint NPAR_P, MSIZE_P, NCH_P, PAR_P, FL_P, I_P, LET_P, CH_P, PAGE_SIZE; 
    nbits = pTrie->nbits;
    setvars(t[1]);
    npar = ((uint)t[NPAR_P])<<1;
    maxsize = N[t[1]].maxsize; 
    if (npar < maxsize) {
       if (npar == 0) i = 1 + PAR_P; 
       else i = ((npar-1)>>3) + 1 + PAR_P;
       if (curbit < npar)
          if (((npar-1)&0x07) == 7) t[i] = (t[i-1] & 0XC0) >>6;
	  else if (((npar-1)&0x07) == 6) t[i] = (t[i-1] & 0X40) >>6;
       for (i--; i > bytepos; i--)
          t[i] = (t[i] << 2) | ((t[i-1] & 0XC0) >> 6);
       /* cubre caso hay insercion al final del bloque y en nuevo byte */
       if ((curbit > 0) && ((curbit & 0x07) == 0) && (curbit >= npar)) i++; 
       t[i] = (((t[i] << 2)
             & (~0 << ((curbit & 0x07) + 2)))
	     | (2 << (curbit & 0x07)))
	     | (t[i] 
	     & ~(~0 << (curbit & 0x07)));
       if ((curbit & 0x07) == 7) t[i+1] |= 1;
       /* Ahora actualizo los flags */
       lastflagbyte = ((cantflags = (npar >> 1))) >> 3;
       curflagbyte  = curflag >> 3; 
       i            = FL_P + lastflagbyte;
       if (curflag < cantflags) { // si la insercion no es al final
          if ((((cantflags-1)&0x07)==7)&&(cantflags > 0)) {
             t[i]=(t[i-1]&0X80)>>7; i--;
          } 
          for (;i > (FL_P + curflagbyte); i--)
             t[i] = (t[i] << 1) | ((t[i-1] & 0X80) >> 7);
	  t[i] = ((t[i] << 1)
                & (~0 << ((curflag & 0x07) + 1)))
                | (t[i]
                & ~(~0 << (curflag & 0x07)));
       }
       /* Ahora inserto la letra en la lista de letras, y el id en su lista */
       i = LET_P + (npar >> 1);
       j = nbits*(npar >> 1);
       if (isrevtrie) {
          for (; i > letterpos; i--) {
             t[i] = t[i-1];
             bitputint((int*)&t[I_P],j,nbits,bitgetint((int*)&t[I_P],
                       j-nbits,nbits));
             j-=nbits;
          } 
	  t[i] = letter;
	  bitputint((int*)&t[I_P],j,nbits,id);
       }
       else {
          for (; i > letterpos; i--) {
             t[i] = t[i-1];
             bitput((uint*)&t[I_P],j,nbits,bitget((uint*)&t[I_P],
                    j-nbits,nbits));
             j-=nbits;
          }
          t[i] = letter;
	  bitput((uint*)&t[I_P],j,nbits,(uint)id);
       } 
       t[NPAR_P]++; // a new node in the page
    }
    else {
       byte **newc, *tnew, *tinsert; int p, curbyte;
       int LET_P_ANT;
       register byte mb, regcond, mbaux;	      
       if (t[1] != lastN) {
          LET_P_ANT = LET_P;
          tinsert   = grow(t, tparent); // crecimiento en la pagina
	  curbyte   = bytepos;
          letterpos = letterpos - LET_P_ANT + N[tinsert[1]].LETTERS_P;
       }
       else { // pagina no puede crecer mas, hay que selecc. y bajar subtrie
          uint cant;
	  tnew      = calloc(N[lastN].Psize, sizeof(byte));
          tnew[1]   = lastN;  // tnew es del maximo tamanio posible  
          p=select_subtrie(t,tnew,curflag,&tinsert,&curbit,pTrie->nbits);
	  //en curbit queda nueva pos de insercion 
          t        = shrink(t, &tinsert);
          *tparent = t;
          setvars(t[1]);
          tnew     = shrink(tnew, &tinsert); // ajusta tnew a su tamanio correcto
          newc     = mallocHeap(pTrie->heaps[t[NCH_P]+1]);
          memcpy(newc, *((byte ***)&t[CH_P]), p*sizeof(byte *));
          memcpy(newc+p+1, (*((byte ***)&t[CH_P]))+p, 
                 (t[NCH_P]-p)*sizeof(byte *));
          if (t[NCH_P]>0) 
             freeHeap(pTrie->heaps[t[NCH_P]], *((byte ***)&t[CH_P]));
          *((byte ***)&t[CH_P])      = newc;
          (*((byte ***)&t[CH_P]))[p] = tnew;
          if (tinsert!=t) {tparent = &((*((byte ***) &t[CH_P]))[p]);}
          t[NCH_P]++; // a new child in t
          setvars(tinsert[1]);
          cant=0;
	  curbyte = PAR_P+(curbit>>3);
          for (i = PAR_P; i < curbyte; i++) cant+= T[(byte)~tinsert[i]];
          cant+= T[(byte)(~tinsert[curbyte]&~(~0<<(curbit&0X07)))];
	  letterpos = LET_P+cant;
	  curflag = cant;
       }       
       /* Llamada recursiva para insertar el nuevo nodo (seguramente el 
          bloque en donde se inserta tiene lugar ==> no entra mas por aca) */
       insert_page(pTrie, tinsert, tparent, curbyte, curbit, letterpos, 
                   curflag, id, letter,isrevtrie);
    }
 }
	

// si el flag indicado por curflag es 1, desciende al hijo
// correpondiente de t

byte *descend(byte *t, uint curflag, bool *success, byte ***tparent,
              uint FL_P, uint CH_P)
 {
    register byte mbaux;
    uint aux, posaux;
    aux = curflag - 1;
    *success = false;
    posaux = (aux >> 3)+FL_P;
    mbaux = 1<<(aux&(uint)0X07);
    if ((t[posaux]&mbaux)!=0) {
       int i, cant = 0;
       // calcula la posicion en el arreglo de hijos a usar en el descenso
       for (i=FL_P; i<posaux; i++) cant+= T[t[i]];
       cant+= T[(t[posaux]&~(~0<<(aux&0X07)))];
       // modifica el puntero al padre
       *tparent = &((*((byte ***)&t[CH_P]))[cant]);
       // desciende al hijo correspondiente
       t = (*((byte ***)&t[CH_P]))[cant];
       *success = true;
    }
    return t;
 }
 
        // inserts word[0...] into pTrie and returns new text ptr
	// insertion proceeds until we get a new trie node

void insertTrie(trie_min pTrie, byte *word, int *m, FILE *fin, FILE *fout)
 { 
    byte *t = pTrie->trie, **tparent=&(pTrie->trie);
    uint bytepos, letterpos, curbit, curflag, npar;
    int excess;
    uint NPAR_P, MSIZE_P, NCH_P, PAR_P, FL_P, I_P, LET_P, CH_P, PAGE_SIZE;
    register byte mb, regcond, mbaux;
    uint i, curid; 
    bool success;
    curflag  = curbit = 0;
    setvars(t[1]);
    npar     = ((uint)t[NPAR_P]) << 1;
    mb       = 1; 
    bytepos  = PAR_P; 
    letterpos = LET_P;
    curid    = 0; //bitget((uint*)&t[I_P],(letterpos-LET_P)*pTrie->nbits,pTrie->nbits);
    if (npar > 0)
    while ((t[bytepos] & mb) == 0) {
       while ((curbit < npar) && ((t[bytepos] & mb) == 0)) {
          if (t[letterpos] >= word[*m]) break;
          excess = 0; mbaux=mb<<1; curbit++; curflag++;
	  while (excess != -1) {
	     regcond = (mbaux==0); mb = mbaux | regcond; bytepos+=regcond;
	     if ((t[bytepos] & mb) == 0) {excess++; curflag++;}
	     else {excess--; letterpos++;}
	     mbaux=mb<<1; curbit++;
	  }
          regcond = (mbaux==0); mb = mbaux|regcond; bytepos+=regcond;
       }
       if ((curbit>=npar)
         || ((t[bytepos]&mb)!=0)
	 || (t[letterpos]>word[*m])) break;
       
       mbaux=mb<<1;
       curbit++; curflag++;
       regcond = (mbaux==0); mb = mbaux|regcond; bytepos+=regcond;
       curid = bitget((uint*)&t[I_P],(letterpos-LET_P)*pTrie->nbits,pTrie->nbits);
       (*m)++; letterpos++;
       if ((*m) == (BUFFER_TEXT_SIZE-1)) {
          txtload(fin, word); 
	  (*m) = 0;
       }
       if ((t[bytepos] & mb)!=0) {//current parenthesis is ')'
          t = descend(t, curflag, &success, &tparent, FL_P, CH_P);
	  if (success) {
	     mb = 1; curbit = 0;
	     curflag = 0;
	     setvars(t[1]);
	     npar = ((uint)t[NPAR_P])<<1;
	     letterpos = LET_P; bytepos = PAR_P;
//             curid = bitget((uint*)&t[I_P],(letterpos-LET_P)*pTrie->nbits,pTrie->nbits);
	  }
       }
    }
    fwrite(&word[*m], sizeof(byte), 1, fout);
    fwrite(&curid, sizeof(uint), 1, fout);   
 
    insert_page(pTrie,t,tparent,bytepos,curbit,letterpos,curflag,
                pTrie->nid++,word[*m],false);
    (*m)++;
 }


static void traverse(lztrie T, byte *t, uint *parent, byte *letters, 
                     uint *pi, uint *pli, uint nbits_semi, 
		     uint nbits_packed, uint *nbytes, int *maxdepth, 
		     int mydepth, FILE *fp) 
// parent: arreglo de parentesis
// letters: arreglo de letras
// ids: arreglo de ids empaquetados
// nbits_semi: # de bits con que se representa cada id en el trie (semi-empaquetado)
// nbits_packed: # de bits con que se representara cada id en el arreglo ids
 {
    register byte mb, mbflag;
    uint npar, curflag, curbytepar, curbyteflag, i, curchild, letterpos, idpos;
    uint NPAR_P, MSIZE_P, NCH_P, PAR_P, FL_P, I_P, LET_P, CH_P, PAGE_SIZE;
    int myid;
    
    setvars(t[1]);
    *nbytes    += PAGE_SIZE + t[NCH_P]*sizeof(byte *);
    npar        = ((uint)t[NPAR_P])<<1;
    mb          = 1; 
    mbflag      = 1;
    curbytepar  = PAR_P; 
    curchild    = 0;   
    curbyteflag = FL_P;
    letterpos   = LET_P; 
    idpos       = 0;
    curflag = 0;
    for (i = 0; i < npar; i++) {
       if ((t[curbytepar]&mb) == 0) {
          mydepth++; 
	  bitclean(parent, *pi);
	  (*pi)++;
          uint id;
          letters[*pli] = t[letterpos];
          id = bitget(((uint *)&t[I_P]),idpos*nbits_semi,nbits_semi);
          fwrite(&id, sizeof(uint), 1, fp);  // writes ids temporarily to disk

	  (*pli)++;
	  letterpos++; idpos++;
	  if ((t[curbyteflag]&mbflag)!=0) {
	     traverse(T,(*((byte ***)&t[CH_P]))[curchild],parent,letters,
	              pi,pli,nbits_semi,nbits_packed,nbytes, maxdepth, mydepth, fp);
	     curchild++;
	  }
	  mbflag <<= 1;
	  curflag++;
	  if (mbflag == 0) {mbflag = 1; curbyteflag++;}
       }
       else {
          mydepth--;
	  bitset(parent, *pi);
	  (*pi)++;
       }
       if (letters) if (mydepth > (*maxdepth)) *maxdepth = mydepth;
       mb <<= 1;
       if (mb == 0) {mb = 1; curbytepar++;}
    }
    free((void *)t);
 }


static void traverseTrie(byte *t, uint *maxdepth, uint mydepth, uint *nbytes)
 {
    register byte mb, mbflag;
    uint npar, curflag, curbytepar, curbyteflag, i, curchild, letterpos, idpos;
    uint NPAR_P, MSIZE_P, NCH_P, PAR_P, FL_P, I_P, LET_P, CH_P, PAGE_SIZE;
    int myid;

    setvars(t[1]);
    *nbytes    += PAGE_SIZE + t[NCH_P]*sizeof(byte *);
    npar        = ((uint)t[NPAR_P])<<1;
    mb          = 1;
    mbflag      = 1;
    curbytepar  = PAR_P;
    curchild    = 0;
    curbyteflag = FL_P;
    letterpos   = LET_P;
    idpos       = 0;
    curflag = 0;
    for (i = 0; i < npar; i++) {
       if ((t[curbytepar]&mb) == 0) {
          mydepth++;
          if ((t[curbyteflag]&mbflag)!=0) {
             traverseTrie((*((byte ***)&t[CH_P]))[curchild], maxdepth, mydepth, nbytes);
             curchild++;
          }
          mbflag <<= 1;
          curflag++;
          if (mbflag == 0) {mbflag = 1; curbyteflag++;}
       }
       else {
          mydepth--;
       }
       if (mydepth > (*maxdepth)) *maxdepth = mydepth;
       mb <<= 1;
       if (mb == 0) {mb = 1; curbytepar++;}
    }
    free((void *)t);
 }





uint representTrie(trie_min pTrie, uint *parent, byte *letters, int *maxdepth, FILE *fp)
 {
    uint pi, pli, nbits;
    uint nbytes = 0;
    nbits = bits(pTrie->nid-1);
    letters[0] = 0; // dummy value
    bitclean(parent,0); // the '(' belonging to the root node
    //bitput(ids,0,nbits,0);
    pi = pli = 1;
    *maxdepth = -1;
    traverse(NULL,pTrie->trie,parent,letters,&pi,&pli,pTrie->nbits,nbits,&nbytes,maxdepth,0,fp);
    bitset(parent,pi); // the ')' belonging to the root node
    printf("Tamanio LZTrie: %d\n",nbytes);
    return nbytes;
 }

uint sizeTrie(trie_min pTrie, uint *maxdepth)
 {  uint nbytes=0;
    traverseTrie(pTrie->trie, maxdepth, 0, &nbytes);
    return nbytes;
 }

	
