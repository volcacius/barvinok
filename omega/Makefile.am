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
	polyfunc.h \
	vertices.cc \
	vertices.h \
	y.tab.cc \
	omega_calc/yylex.h
if USE_PARKER
    PARKER_LA = ../parker/libparker.la
endif
LDADD = $(PARKER_LA) \
	-lcode_gen -lomega ../libbarvinok.la ../bernstein/libbernstein.la \
	@POLYLIB_LIBS@ @PIPLIB_LIBS@
AM_LDFLAGS = @OMEGA_LDFLAGS@ @POLYLIB_LDFLAGS@ @PIPLIB_LDFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/bernstein/include \
	    @OMEGA_CPPFLAGS@

EXTRA_DIST = \
	parser.l \
	parser.y

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<

