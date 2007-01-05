if HAVE_GINAC
    MAYBE_BERNSTEIN = bernstein
endif
SUBDIRS = $(MAYBE_BERNSTEIN) lib . omega doc

AM_CPPFLAGS = -I$(srcdir)/bernstein/include -I$(srcdir)/lib
AM_LDFLAGS = -Llib

lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count randomtest barvinok_enumerate \
		  barvinok_ehrhart \
		  verif_ehrhart barvinok_enumerate_e \
		  barvinok_series remove_redundant_equalities \
		  barvinok_union \
		  @bv_extra_programs@
EXTRA_PROGRAMS = piptest verify_lexsmaller polyhedron_sample 4coins lexmin \
		 polytope_scan
pkginclude_HEADERS = \
    barvinok/NTL_QQ.h \
    barvinok/barvinok.h \
    barvinok/util.h \
    barvinok/evalue.h \
    barvinok/genfun.h \
    barvinok/options.h \
    barvinok/polylib.h \
    barvinok/sample.h \
    barvinok/basis_reduction.h \
    barvinok/bernstein.h
LINK = $(CXXLINK)

if HAVE_GINAC
    BERNSTEIN_CC = bernstein.cc
endif
if HAVE_CDDLIB
    BR_CDD = basis_reduction_cdd.cc initcdd.cc
endif
if HAVE_GLPK
    BR_GLPK = basis_reduction_glpk.c
endif
libbarvinok_la_SOURCES = \
    barvinok/evalue.h \
    barvinok/genfun.h \
    barvinok/util.h \
    barvinok/barvinok.h \
    NTL_QQ.cc \
    basis_reduction.c \
    evalue.c \
    genfun.cc \
    util.c \
    version.c \
    bfcounter.cc \
    bfcounter.h \
    conversion.cc \
    conversion.h \
    decomposer.cc \
    decomposer.h \
    dpoly.cc \
    dpoly.h \
    genfun_constructor.cc \
    genfun_constructor.h \
    lattice_point.cc \
    lattice_point.h \
    options.c \
    reduce_domain.c \
    reduce_domain.h \
    reducer.cc \
    reducer.h \
    remove_equalities.h \
    scarf.cc \
    mat_util.cc \
    mat_util.h \
    barvinok.cc \
    $(BERNSTEIN_CC) \
    $(BR_CDD) \
    $(BR_GLPK)
EXTRA_libbarvinok_la_SOURCES = \
    bernstein.cc \
    piputil.h \
    initcdd.cc \
    initcdd.h
if HAVE_GINAC
    BERNSTEIN_LA = $(top_builddir)/bernstein/libbernstein.la
endif
libbarvinok_la_LIBADD = @LTLIBOBJS@ $(BERNSTEIN_LA) lib/libgnu.la
libbarvinok_la_LDFLAGS = @BV_LDFLAGS@ -version-info @versioninfo@
LDADD = libbarvinok.la

test_SOURCES = test.c
barvinok_count_SOURCES = barvinok_count.c
barvinok_ehrhart_SOURCES = barvinok_ehrhart.cc
barvinok_series_SOURCES = barvinok_series.cc
barvinok_union_SOURCES = barvinok_union.cc
if HAVE_OMEGA
BEEO_SOURCES = omega/Exit.cc omega/convert.cc
else
BEEO_SOURCES = 
endif
barvinok_enumerate_e_SOURCES = \
	verify.h \
	verify.c \
	barvinok_enumerate_e.cc \
	$(BEEO_SOURCES)
barvinok_enumerate_e_CPPFLAGS = @OMEGA_CPPFLAGS@ $(AM_CPPFLAGS)
barvinok_enumerate_e_LDFLAGS = @OMEGA_LDFLAGS@ @LDFLAGS@
barvinok_enumerate_e_LDADD = @OMEGA_LIBS@ libbarvinok.la
verif_ehrhart_SOURCES = verif_ehrhart.c verify_main.cc verif_ehrhart.h
4coins_SOURCES = 4coins.cc
if HAVE_CDDLIB
POLYSIGN = polysign_cddf.cc polysign_cdd.cc
endif
lexmin_SOURCES = \
	lexmin.h \
	lexmin.cc \
	combine.c \
	combine.h \
	edomain.cc \
	edomain.h \
	evalue_util.cc \
	evalue_util.h \
	fdstream.cc \
	fdstream.h \
	$(POLYSIGN) \
	polysign.c \
	polysign_polylib.c \
	polysign.h \
	verify.h \
	verify.c
EXTRA_lexmin_SOURCE = polysign_cddf.cc polysign_cdd.cc polysign.c

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = \
    ChangeLog \
    $(TESTFILES) \
    latte2polylib.pl \
    NTL_5_3_2.patch \
    doc/Internal.tex \
    doc/applications.tex \
    doc/barvinok.tex \
    doc/omega.tex \
    doc/polymake.tex \
    doc/reports.tex \
    doc/barvinok.bib \
    doc/barvinok.gdf \
    doc/chicago.bst \
    doc/chicago.sty \
    doc/mydefs.sty \
    polymake/Makefile.in \
    polymake/README \
    polymake/autogen.sh \
    polymake/configure.in \
    polymake/lattice_points.cc

ACLOCAL_AMFLAGS = -I m4

install-data-local:
	@test -z "$(pkgconfig_libdir)" || $(mkdir_p) "$(DESTDIR)$(pkgconfig_libdir)"
	$(INSTALL_DATA) $(pkgconfig_libfile) "$(DESTDIR)$(pkgconfig_libdir)/$(pkgconfig_libfile)"

uninstall-local:
	rm -f "$(DESTDIR)$(pkgconfig_libdir)/$(pkgconfig_libfile)"

dist-hook:
	(cd $(distdir)/polymake; \
	    ./autogen.sh; rm -rf autogen.sh autom4te.cache)
	(cd doc; make barvinok.pdf)
	cp doc/barvinok.pdf $(distdir)/doc/
	if test -f $(top_srcdir)/.git/HEAD; then \
	    echo @GIT_HEAD_VERSION@ > $(distdir)/GIT_HEAD_ID; \
	else \
	    echo $(GIT_HEAD_ID) > $(distdir)/GIT_HEAD_ID; \
	fi

if HAVE_PIPLIB
BEE_TESTDIRS = ehrhart_e ehrhart_e/scarf ehrhart_e/piplib
else
BEE_TESTDIRS = ehrhart_e ehrhart_e/scarf
endif

check-series:
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i; \
		./verif_ehrhart --series < $$i; \
	    fi \
	done
check:
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i; \
		./verif_ehrhart < $$i; \
	    fi \
	done
	@for dir in $(BEE_TESTDIRS); do \
	    for options in '' '--pip' '--pip --omega'; do \
		for i in $(top_srcdir)/tests/$$dir/*; do \
		    if test -f $$i; then \
			echo $$i $$options; \
			./barvinok_enumerate_e --verify $$options < $$i; \
		    fi \
		done \
	    done \
	done
	@for i in $(top_srcdir)/tests/ehrhart_e/scarf/*; do \
	    if test -f $$i; then \
		echo $$i --scarf; \
		./barvinok_enumerate_e --verify --scarf < $$i; \
	    fi \
	done
	@for i in $(top_srcdir)/tests/*; do \
	    for options in '' '--primal'; do \
		if test -f $$i; then \
		    echo -n $$i $$options; \
		    ./test $$options < $$i || echo -n " NOT"; \
		    echo " ok"; \
		fi \
	    done \
	done
if HAVE_GLPK
	@for i in $(top_srcdir)/tests/lexmin/*; do \
	    if test -f $$i; then \
		echo $$i; \
		./lexmin --verify < $$i; \
	    fi \
	done
endif

version.h: @GIT_HEAD@
	echo '#define GIT_HEAD_ID "'@GIT_HEAD_VERSION@'"' > $@
