MONAPATH = /home/skimo/src/mona-1.4

CC = gcc
CFLAGS = -O3 -DNDEBUG

MEM_OBJ = $(MONAPATH)/Mem/libmem.a
BDD_OBJ = $(MONAPATH)/BDD/libbdd.a
DFA_OBJ = $(MONAPATH)/DFA/libdfa.a

INCLUDES = -I$(MONAPATH)/Mem -I$(MONAPATH)/DFA -I$(MONAPATH)/BDD 

count: 
	$(CC) $(CFLAGS) $(INCLUDES) -o count $(in).c construction.c count_paths.c $(DFA_OBJ) $(BDD_OBJ) $(MEM_OBJ) -lm
