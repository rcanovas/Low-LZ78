CC = gcc 
FLAGS = -m32 -O9 -lm 

##############################
#  Objets  (.o)             #
##############################

OBJ_UnCom = trie_uc.o heap.o basics.o

OBJ_Min = lztrie.o trie_min.o parentheses.o bitmap.o hash.o heap.o basics.o

##############################
# make all                   #
##############################

all: all_obj clz78_uc ulz78_uc clz78_min ulz78_min


all_obj:  lztrie trie_min trie_uc parentheses bitmap hash heap basics

###############################
# crear objetos y ejecutables #
###############################

clz78_uc:
	$(CC) $(FLAGS) $(OBJ_UnCom) -o Clz78_uc clz78_uc.c

clz78_min:
	$(CC) $(FLAGS) $(OBJ_Min) -o Clz78_min clz78_min.c

ulz78_uc:
	$(CC) $(FLAGS) $(OBJ_UnCom) -o Dlz78_uc ulz78_uc.c

ulz78_min:
	$(CC) $(FLAGS) $(OBJ_UnCom) -o Dlz78_min ulz78_min.c

trie_uc:  
	$(CC) $(FLAGS) -c -o trie_uc.o trie_uc.c

trie_min:  
	$(CC) $(FLAGS) -c -o trie_min.o trie_min.c

heap:  	
	$(CC) $(FLAGS) -c -o heap.o heap.c

hash:  	
	$(CC) $(FLAGS) -c -o hash.o hash.c

basics:
	$(CC) $(FLAGS) -c -o basics.o basics.c

bitmap:
	$(CC) $(FLAGS) -c -o bitmap.o bitmap.c

parentheses:
	$(CC) $(FLAGS) -c -o parentheses.o parentheses.c

lztrie:  
	$(CC) $(FLAGS) -c -o lztrie.o lztrie.c

##############################
# clean                      #
##############################
clean:
	  rm -f *.o *.c~ *.h~ makefile~ Clz78_* Dlz78_*
