if HAVE_GINAC
    MAYBE_BERNSTEIN = bernstein
endif
if BUNDLED_POLYLIB
    MAYBE_POLYLIB = polylib
endif
if BUNDLED_PIPLIB
    MAYBE_PIPLIB = piplib
endif
SUBDIRS = $(MAYBE_PIPLIB) $(MAYBE_POLYLIB) $(MAYBE_BERNSTEIN) lib . omega doc

FORCE:
lib/libgnu.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C lib libgnu.la
bernstein/libbernstein.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C bernstein libbernstein.la
polylib/libpolylibgmp.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C polylib libpolylibgmp.la
piplib/libpiplibMP.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C piplib libpiplibMP.la

AM_CPPFLAGS = -I$(srcdir)/bernstein/include -I$(srcdir)/lib \
	      @POLYLIB_CPPFLAGS@ @PIPLIB_CPPFLAGS@ @GINACLIB_CPPFLAGS@
AM_LDFLAGS = -Llib @POLYLIB_LDFLAGS@ @PIPLIB_LDFLAGS@

lib_LTLIBRARIES = libbarvinok.la
bin_PROGRAMS = barvinok_count barvinok_enumerate barvinok_enumerate_e \
		  barvinok_ehrhart polytope_scan
noinst_PROGRAMS = test testlib randomtest \
		  remove_redundant_equalities \
		  barvinok_union polytope_volume test_approx \
		  barvinok_summate verify_lexsmaller piptest \
		  polyhedron_sample 4coins lexmin @bv_barvinok_maximize@
EXTRA_PROGRAMS = barvinok_maximize
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
    POLYSIGN_CDD = polysign_cddf.cc polysign_cdd.cc
endif
if HAVE_GLPK
    BR_GLPK = basis_reduction_glpk.c
    POLYSIGN_GLPK = polysign_glpk.c
endif
if HAVE_TOPCOM
    TOPCOM = topcom.c topcom.h
endif
libbarvinok_la_SOURCES = \
    barvinok/evalue.h \
    barvinok/genfun.h \
    barvinok/util.h \
    barvinok/barvinok.h \
    NTL_QQ.cc \
    basis_reduction.c \
    basis_reduction_pip.c \
    basis_reduction_pip_dual.c \
    evalue.c \
    genfun.cc \
    util.c \
    version.c \
    bernoulli.c \
    bernoulli.h \
    bfcounter.cc \
    bfcounter.h \
    conversion.cc \
    conversion.h \
    counter.cc \
    counter.h \
    decomposer.cc \
    decomposer.h \
    dpoly.cc \
    dpoly.h \
    euler.cc \
    euler.h \
    genfun_constructor.cc \
    genfun_constructor.h \
    lattice_point.cc \
    lattice_point.h \
    options.c \
    param_util.c \
    param_util.h \
    piputil.h \
    piputil.c \
    $(POLYSIGN_CDD) \
    $(POLYSIGN_GLPK) \
    polysign.c \
    polysign_polylib.c \
    polysign.h \
    reduce_domain.c \
    reduce_domain.h \
    reducer.cc \
    reducer.h \
    remove_equalities.h \
    remove_equalities.c \
    sample.c \
    scale.c \
    scale.h \
    scarf.cc \
    series.cc \
    tcounter.cc \
    tcounter.h \
    $(TOPCOM) \
    volume.c \
    volume.h \
    mat_util.cc \
    mat_util.h \
    matrix_read.cc \
    matrix_read.h \
    barvinok.cc \
    barvinok_e.cc \
    $(BERNSTEIN_CC) \
    $(BR_CDD) \
    $(BR_GLPK)
EXTRA_libbarvinok_la_SOURCES = \
    basis_reduction_cdd.cc \
    basis_reduction_glpk.c \
    bernstein.cc \
    polysign_cddf.cc \
    polysign_cdd.cc \
    polysign.c \
    initcdd.cc \
    initcdd.h \
    topcom.c \
    topcom.h
if HAVE_GINAC
    BERNSTEIN_LA = $(top_builddir)/bernstein/libbernstein.la
endif
if BUNDLED_POLYLIB
    POLYLIB_LA = $(top_builddir)/polylib/libpolylibgmp.la
endif
if BUNDLED_PIPLIB
    PIPLIB_LA = $(top_builddir)/piplib/libpiplibMP.la
endif
libbarvinok_la_LIBADD = @LTLIBOBJS@ $(POLYLIB_LA) @POLYLIB_LIBS@ \
			$(PIPLIB_LA) @PIPLIB_LIBS@ $(BERNSTEIN_LA) \
			lib/libgnu.la
libbarvinok_la_LDFLAGS = @BV_LDFLAGS@ -version-info @versioninfo@
LDADD = libbarvinok.la @POLYLIB_LIBS@ @PIPLIB_LIBS@

test_SOURCES = test.c
testlib_SOURCES = \
	evalue_read.c \
	evalue_read.h \
	testlib.cc
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
	skewed_genfun.cc \
	skewed_genfun.h \
	verify.h \
	verify.c \
	verif_ehrhart.h \
	verif_ehrhart.c \
	verify_series.h \
	verify_series.cc
EXTRA_barvinok_enumerate_SOURCES = fdstream.cc fdstream.h
barvinok_enumerate_e_SOURCES = \
	barvinok_enumerate_e.cc \
	evalue_convert.cc \
	evalue_convert.h \
	$(FDSTREAM) \
	skewed_genfun.cc \
	skewed_genfun.h \
	verify.h \
	verify.c \
	verif_ehrhart.h \
	verif_ehrhart.c \
	verify_series.h \
	verify_series.cc \
	$(BEEO_SOURCES)
EXTRA_barvinok_enumerate_e_SOURCES = fdstream.cc fdstream.h
barvinok_enumerate_e_CPPFLAGS = @OMEGA_CPPFLAGS@ $(AM_CPPFLAGS)
barvinok_enumerate_e_LDFLAGS = @OMEGA_LDFLAGS@ @LDFLAGS@
barvinok_enumerate_e_LDADD = @OMEGA_LIBS@ libbarvinok.la @POLYLIB_LIBS@ @PIPLIB_LIBS@
barvinok_maximize_SOURCES = \
	evalue_convert.cc \
	evalue_convert.h \
	evalue_read.c \
	evalue_read.h \
	$(FDSTREAM) \
	maximize.cc \
	verify.h \
	verify.c
EXTRA_barvinok_maximize_SOURCES = fdstream.cc fdstream.h
barvinok_summate_SOURCES = \
	evalue_convert.cc \
	evalue_convert.h \
	evalue_read.c \
	evalue_read.h \
	$(FDSTREAM) \
	summate.cc \
	verify.h \
	verify.c
EXTRA_barvinok_summate_SOURCES = fdstream.cc fdstream.h
4coins_SOURCES = 4coins.cc
lexmin_SOURCES = \
	lexmin.h \
	lexmin.cc \
	combine.c \
	combine.h \
	edomain.cc \
	edomain.h \
	evalue_util.cc \
	evalue_util.h \
	verify.h \
	verify.c
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
    polymake/convert.cc \
    polymake/convert.h \
    polymake/h_star_vector.cc \
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

dist-git: dist
	mv $(distdir).tar.gz @GIT_HEAD_VERSION@.tar.gz

BEE_TESTDIRS = ehrhart_e ehrhart_e/scarf ehrhart_e/piplib

check: check-testlib check-enumerate check-enumerate_e check-test \
	check-lexmin check-approx check-evalue check-euler
check-testlib: testlib$(EXEEXT)
	./testlib$(EXEEXT)
check-approx: test_approx$(EXEEXT)
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i | ./test_approx$(EXEEXT) -q -r 10 || exit; \
	    fi \
	done
check-evalue: @bv_barvinok_maximize@ barvinok_summate$(EXEEXT)
	@for i in $(top_srcdir)/tests/evalue/*; do \
	    if test -f $$i; then \
		echo $$i; \
		if test -n "@bv_barvinok_maximize@"; then \
		    ./barvinok_maximize$(EXEEXT) -T < $$i || exit; \
		fi; \
		./barvinok_summate$(EXEEXT) -T < $$i || exit; \
	    fi \
	done
check-euler: barvinok_summate$(EXEEXT)
	@for i in $(top_srcdir)/tests/euler/*; do \
	    if test -f $$i; then \
		echo $$i; \
		./barvinok_summate$(EXEEXT) -T --summation=euler < $$i || exit; \
	    fi \
	done
if HAVE_TOPCOM
    TOPCOM_CD = '--chamber-decomposition=topcom'
endif
check-enumerate: barvinok_enumerate$(EXEEXT)
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i; \
		for options in '' '--series' '--series --primal' $(TOPCOM_CD); do \
		    for spec in 'random' 'bf' 'df'; do \
			opt="--specialization=$$spec $$options"; \
			echo "        $$opt"; \
			./barvinok_enumerate$(EXEEXT) --verify $$opt < $$i || exit; \
		    done \
		done \
	    fi \
	done
check-enumerate_e: barvinok_enumerate_e$(EXEEXT)
	@for dir in $(BEE_TESTDIRS); do \
	    for i in $(top_srcdir)/tests/$$dir/*; do \
		if test -f $$i; then \
		    for options in '' '--pip' '--pip --omega'; do \
			for spec in 'random' 'bf' 'df'; do \
			    opt="--specialization=$$spec $$options"; \
			    echo $$i $$opt; \
			    ./barvinok_enumerate_e$(EXEEXT) --verify $$opt < $$i || exit; \
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
		    ./barvinok_enumerate_e$(EXEEXT) --verify $$opt < $$i || exit; \
		done \
	    fi \
	done
check-test: test$(EXEEXT)
	@for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		for options in '--index=10' '--primal --index=10'; do \
		    for spec in 'random' 'bf' 'df' 'todd'; do \
			opt="--specialization=$$spec $$options"; \
			echo -n $$i $$opt; \
			./test$(EXEEXT) $$opt < $$i || echo -n " NOT"; \
			echo " ok"; \
		    done \
		done \
	    fi \
	done
if HAVE_GLPK
check-lexmin: lexmin$(EXEEXT)
	@for i in $(top_srcdir)/tests/lexmin/*; do \
	    if test -f $$i; then \
		for spec in 'random' 'bf' 'df'; do \
		    opt="--specialization=$$spec"; \
		    echo $$i $$opt; \
		    ./lexmin$(EXEEXT) --verify $$opt < $$i || exit; \
		done \
	    fi \
	done
else
check-lexmin:
endif

version.h: @GIT_HEAD@
	echo '#define GIT_HEAD_ID "'@GIT_HEAD_VERSION@'"' > $@
