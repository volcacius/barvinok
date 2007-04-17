if HAVE_GINAC
    MAYBE_BERNSTEIN = bernstein
endif
SUBDIRS = $(MAYBE_BERNSTEIN) lib . omega doc

lib/libgnu.la:
	$(MAKE) $(AM_MAKEFLAGS) -C lib libgnu.la
bernstein/libbernstein.la:
	$(MAKE) $(AM_MAKEFLAGS) -C bernstein libbernstein.la

AM_CPPFLAGS = -I$(srcdir)/bernstein/include -I$(srcdir)/lib @GINACLIB_CPPFLAGS@
AM_LDFLAGS = -Llib

lib_LTLIBRARIES = libbarvinok.la
bin_PROGRAMS = barvinok_count barvinok_enumerate barvinok_enumerate_e \
		  barvinok_ehrhart
noinst_PROGRAMS = test randomtest \
		  remove_redundant_equalities \
		  barvinok_union polytope_volume test_approx \
		  @bv_extra_programs@
EXTRA_PROGRAMS = piptest verify_lexsmaller polyhedron_sample 4coins lexmin \
		 polytope_scan barvinok_maximize
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
nodist_pkginclude_HEADERS = \
    barvinok/NTL.h \
    barvinok/set.h
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
    param_util.c \
    param_util.h \
    reduce_domain.c \
    reduce_domain.h \
    reducer.cc \
    reducer.h \
    remove_equalities.h \
    scale.c \
    scale.h \
    scarf.cc \
    volume.c \
    volume.h \
    mat_util.cc \
    mat_util.h \
    barvinok.cc \
    $(BERNSTEIN_CC) \
    $(BR_CDD) \
    $(BR_GLPK)
EXTRA_libbarvinok_la_SOURCES = \
    basis_reduction_cdd.cc \
    basis_reduction_glpk.c \
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
barvinok_union_SOURCES = barvinok_union.cc
if HAVE_OMEGA
BEEO_SOURCES = omega/Exit.cc omega/convert.cc
else
BEEO_SOURCES = 
endif
if USE_FDSTREAM
FDSTREAM = fdstream.cc fdstream.h
endif
barvinok_enumerate_SOURCES = \
	barvinok_enumerate.cc \
	evalue_convert.cc \
	evalue_convert.h \
	$(FDSTREAM) \
	verify.h \
	verify.c \
	verif_ehrhart.h \
	verif_ehrhart.c
EXTRA_barvinok_enumerate_SOURCES = fdstream.cc fdstream.h
barvinok_enumerate_e_SOURCES = \
	barvinok_enumerate_e.cc \
	evalue_convert.cc \
	evalue_convert.h \
	$(FDSTREAM) \
	verify.h \
	verify.c \
	verif_ehrhart.h \
	verif_ehrhart.c \
	$(BEEO_SOURCES)
EXTRA_barvinok_enumerate_e_SOURCES = fdstream.cc fdstream.h
barvinok_enumerate_e_CPPFLAGS = @OMEGA_CPPFLAGS@ $(AM_CPPFLAGS)
barvinok_enumerate_e_LDFLAGS = @OMEGA_LDFLAGS@ @LDFLAGS@
barvinok_enumerate_e_LDADD = @OMEGA_LIBS@ libbarvinok.la
barvinok_maximize_SOURCES = \
	evalue_convert.cc \
	evalue_convert.h \
	$(FDSTREAM) \
	maximize.cc \
	verify.h \
	verify.c
EXTRA_barvinok_maximize_SOURCES = fdstream.cc fdstream.h
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
	$(POLYSIGN) \
	polysign.c \
	polysign_polylib.c \
	polysign.h \
	verify.h \
	verify.c
EXTRA_lexmin_SOURCES = polysign_cddf.cc polysign_cdd.cc polysign.c
test_approx_SOURCES = \
	test_approx.c \
	verify.h \
	verify.c

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = \
    ChangeLog \
    $(TESTFILES) \
    latte2polylib.pl \
    NTL_5_3_2.patch \
    basis_reduction_templ.c \
    polysign_cdd_template.cc \
    barvinok/NTL.h.broken \
    barvinok/NTL.h.normal \
    barvinok/set.h.broken \
    barvinok/set.h.normal \
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

check: check-enumerate check-enumerate_e check-test check-lexmin check-approx
check-approx: test_approx
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i | ./test_approx -q -r 10 ; \
	    fi \
	done
check-enumerate: barvinok_enumerate
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i; \
		for options in '' '--series' '--series --primal'; do \
		    for spec in 'random' 'bf' 'df'; do \
			opt="--specialization=$$spec $$options"; \
			echo "        $$opt"; \
			./barvinok_enumerate --verify $$opt < $$i; \
		    done \
		done \
	    fi \
	done
check-enumerate_e: barvinok_enumerate_e
	@for dir in $(BEE_TESTDIRS); do \
	    for i in $(top_srcdir)/tests/$$dir/*; do \
		if test -f $$i; then \
		    for options in '' '--pip' '--pip --omega'; do \
			for spec in 'random' 'bf' 'df'; do \
			    opt="--specialization=$$spec $$options"; \
			    echo $$i $$opt; \
			    ./barvinok_enumerate_e --verify $$opt < $$i; \
			done \
		    done \
		fi \
	    done \
	done
	@for i in $(top_srcdir)/tests/ehrhart_e/scarf/*; do \
	    if test -f $$i; then \
		for spec in 'random' 'bf' 'df'; do \
		    opt="--specialization=$$spec --scarf"; \
		    echo $$i $$opt; \
		    ./barvinok_enumerate_e --verify $$opt < $$i; \
		done \
	    fi \
	done
check-test: test
	@for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		for options in '--index=10' '--primal --index=10'; do \
		    for spec in 'random' 'bf' 'df'; do \
			opt="--specialization=$$spec $$options"; \
			echo -n $$i $$opt; \
			./test $$opt < $$i || echo -n " NOT"; \
			echo " ok"; \
		    done \
		done \
	    fi \
	done
if HAVE_GLPK
check-lexmin: lexmin
	@for i in $(top_srcdir)/tests/lexmin/*; do \
	    if test -f $$i; then \
		for spec in 'random' 'bf' 'df'; do \
		    opt="--specialization=$$spec"; \
		    echo $$i $$opt; \
		    ./lexmin --verify $$opt < $$i; \
		done \
	    fi \
	done
else
check-lexmin:
endif

version.h: @GIT_HEAD@
	echo '#define GIT_HEAD_ID "'@GIT_HEAD_VERSION@'"' > $@
