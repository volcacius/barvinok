noinst_PROGRAMS = @bv_omega_programs@
EXTRA_PROGRAMS = occ

if BUNDLED_OMEGA
    OMEGA_CODEGEN_LA = $(top_builddir)/omega/omega_lib/libomega.la \
		       $(top_builddir)/omega/code_gen/libcode_gen.la
endif

FORCE:
$(top_builddir)/omega/omega_lib/libomega.la: FORCE
	cd $(top_builddir)/omega/omega_lib; $(MAKE) $(AM_MAKEFLAGS) libomega.la
$(top_builddir)/omega/code_gen/libcode_gen.la: FORCE
	cd $(top_builddir)/omega/code_gen; $(MAKE) $(AM_MAKEFLAGS) libcode_gen.la

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
	@OMEGA_CODEGEN_LIBS@ $(OMEGA_CODEGEN_LA) ../libbarvinok.la
AM_LDFLAGS = @OMEGA_LDFLAGS@ @POLYLIB_LDFLAGS@ @PIPLIB_LDFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/bernstein/include \
	    @POLYLIB_CPPFLAGS@ @OMEGA_CPPFLAGS@ @NTL_CPPFLAGS@

EXTRA_DIST = \
	parser.l \
	parser.y

YACC_FLAGS = -d

y.tab.h y.tab.cc: parser.y
	$(YACC) -v -d $<
	mv y.tab.c y.tab.cc

lex.yy.cc: parser.l y.tab.h
	$(LEX) -olex.yy.cc -i $<

