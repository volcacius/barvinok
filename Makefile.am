MONAPATH = @MONAPATH@
MONA_OBJ = @MONA_OBJ@

lib_LTLIBRARIES = libparker.la
libparker_la_SOURCES = \
	construction.c \
	count_paths.c \
	count_solutions.cc
libparker_la_LIBADD = $(MONA_OBJ)

noinst_PROGRAMS = occ
occ_SOURCES = \
	Exit.cc \
	lex.yy.cc \
	y.tab.cc
LDADD = -lcode_gen -lomega libparker.la

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<

CFLAGS = -O3 -DNDEBUG

INCLUDES = -I$(MONAPATH)/Mem -I$(MONAPATH)/DFA -I$(MONAPATH)/BDD 

AM_CPPFLAGS = $(INCLUDES)
