MONAPATH = @MONAPATH@
MONA_OBJ = @MONA_OBJ@

noinst_PROGRAMS = occ
occ_SOURCES = \
	Exit.cc \
	lex.yy.cc \
	y.tab.cc \
	construction.c \
	count_paths.c \
	count_solutions.cc
LDADD = -lcode_gen -lomega $(MONA_OBJ)

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<

CFLAGS = -O3 -DNDEBUG

INCLUDES = -I$(MONAPATH)/Mem -I$(MONAPATH)/DFA -I$(MONAPATH)/BDD 

AM_CPPFLAGS = $(INCLUDES)
