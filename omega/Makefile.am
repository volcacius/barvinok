noinst_PROGRAMS = @bv_omega_programs@
EXTRA_PROGRAMS = occ

occ_SOURCES = \
	Exit.cc \
	convert.cc \
	convert.h \
	count.cc \
	count.h \
	lex.yy.cc \
	polyfunc.cc \
	y.tab.cc \
	omega_calc/yylex.h
LDADD = -lcode_gen -lomega ../libbarvinok.la ../bernstein/libbernstein.la
LDFLAGS = @OMEGA_LDFLAGS@ @LDFLAGS@
CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/bernstein/include \
	    @OMEGA_CPPFLAGS@ @CPPFLAGS@

EXTRA_DIST = \
	parser.l \
	parser.y

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<

