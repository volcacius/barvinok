noinst_PROGRAMS = occ
occ_SOURCES = \
	Exit.cc \
	lex.yy.cc \
	y.tab.cc \
	construction.c \
	count_paths.c \
	count_solutions.cc
LDADD = -lcode_gen -lomega $(DFA_OBJ) $(BDD_OBJ) $(MEM_OBJ)

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<

MONAPATH = /home/skimo/src/mona-1.4

CC = gcc
CFLAGS = -O3 -DNDEBUG

MEM_OBJ = $(MONAPATH)/Mem/libmem.a
BDD_OBJ = $(MONAPATH)/BDD/libbdd.a
DFA_OBJ = $(MONAPATH)/DFA/libdfa.a

INCLUDES = -I$(MONAPATH)/Mem -I$(MONAPATH)/DFA -I$(MONAPATH)/BDD 

AM_CPPFLAGS = $(INCLUDES)

count: count_solutions.o Exit.o
	$(CC) $(CFLAGS) $(INCLUDES) -o count Exit.o count_solutions.o construction.c count_paths.c $(DFA_OBJ) $(BDD_OBJ) $(MEM_OBJ) -lm -lomega -lstdc++
