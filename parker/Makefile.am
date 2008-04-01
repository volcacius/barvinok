MONAPATH = @MONAPATH@

lib_LTLIBRARIES = libparker.la
libparker_la_SOURCES = \
	construction.c \
	count_paths.c \
	count_solutions.cc \
	count_solutions.h
include_HEADERS = count_solutions.h
libparker_la_LIBADD = -ldfa -lbdd -lmem
libparker_la_LDFLAGS = -L$(MONAPATH)/DFA -L$(MONAPATH)/BDD -L$(MONAPATH)/Mem
INCLUDES = -I$(MONAPATH)/Mem -I$(MONAPATH)/DFA -I$(MONAPATH)/BDD
libparker_la_CPPFLAGS = $(INCLUDES) @OMEGA_CPPFLAGS@ -DNDEBUG

noinst_PROGRAMS = occ
occ_SOURCES = \
	Exit.cc \
	lex.yy.cc \
	y.tab.cc
occ_CPPFLAGS = @OMEGA_CPPFLAGS@
occ_LDFLAGS = @OMEGA_LDFLAGS@
occ_LDADD = -lcode_gen -lomega libparker.la

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<
