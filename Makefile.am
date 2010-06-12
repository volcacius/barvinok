if BUNDLED_POLYLIB
    MAYBE_POLYLIB = polylib
endif
if BUNDLED_PIPLIB
    MAYBE_PIPLIB = piplib
endif
if BUNDLED_ISL
    MAYBE_ISL = isl
    ISL_LA = $(top_builddir)/isl/libisl.la
endif
if BUNDLED_ISL_POLYLIB
    MAYBE_ISL_POLYLIB = isl-polylib
    ISL_POLYLIB_LA = $(top_builddir)/isl-polylib/libisl-polylib.la
endif
if BUNDLED_OMEGA
    MAYBE_OMEGA = omega
    OMEGA_LA = $(top_builddir)/omega/omega_lib/libomega.la
endif
if USE_ZSOLVE
    MAYBE_ZSOLVE = zsolve
endif
if USE_PARKER
    MAYBE_PARKER = parker
endif
SUBDIRS = $(MAYBE_PIPLIB) $(MAYBE_POLYLIB) $(MAYBE_ISL) $(MAYBE_ISL_POLYLIB) \
		$(MAYBE_OMEGA) \
		$(MAYBE_ZSOLVE) $(MAYBE_PARKER) . doc
DIST_SUBDIRS = $(MAYBE_PIPLIB) $(MAYBE_POLYLIB) $(MAYBE_ISL) \
		$(MAYBE_ISL_POLYLIB) $(MAYBE_OMEGA) \
		zsolve parker doc

FORCE:
polylib/libpolylibgmp.la: FORCE
	cd polylib; $(MAKE) $(AM_MAKEFLAGS) libpolylibgmp.la
piplib/libpiplibMP.la: FORCE
	cd piplib; $(MAKE) $(AM_MAKEFLAGS) libpiplibMP.la
isl/libisl.la: FORCE
	cd isl; $(MAKE) $(AM_MAKEFLAGS) libisl.la
isl-polylib/libisl-polylib.la: FORCE
	cd isl-polylib; $(MAKE) $(AM_MAKEFLAGS) libisl-polylib.la
omega/omega_lib/libomega.la: FORCE
	cd omega/omega_lib; $(MAKE) $(AM_MAKEFLAGS) libomega.la
zsolve/libzsolve.la: FORCE
	cd zsolve; $(MAKE) $(AM_MAKEFLAGS) libzsolve.la
parker/libparker.la: FORCE
	cd parker; $(MAKE) $(AM_MAKEFLAGS) libparker.la

AM_CPPFLAGS = @POLYLIB_CPPFLAGS@ @PIPLIB_CPPFLAGS@ \
	      @NTL_CPPFLAGS@ @ISL_CFLAGS@ @ISL_POLYLIB_CFLAGS@
AM_LDFLAGS = @POLYLIB_LDFLAGS@ @PIPLIB_LDFLAGS@ @NTL_LDFLAGS@

LIB_ISL = $(ISL_LA) @ISL_LIBS@
LIB_ISL_POLYLIB = $(ISL_POLYLIB_LA) @ISL_POLYLIB_LIBS@

lib_LTLIBRARIES = libbarvinok.la
bin_PROGRAMS = barvinok_count barvinok_enumerate barvinok_enumerate_e \
		  barvinok_ehrhart polytope_scan iscc
noinst_PROGRAMS = test testlib randomtest \
		  remove_redundant_equalities \
		  barvinok_union polytope_volume test_approx \
		  barvinok_summate verify_lexsmaller \
		  polyhedron_sample 4coins lexmin \
		  @bv_cone_hilbert_basis@ cone_integer_hull \
		  polytope_lattice_width polytope_minimize \
		  polyhedron_integer_hull vector_partition_chambers \
		  semigroup_holes evalue_convert \
		  barvinok_bound test_bound
EXTRA_PROGRAMS = cone_hilbert_basis
pkginclude_HEADERS = \
    barvinok/NTL_QQ.h \
    barvinok/barvinok.h \
    barvinok/util.h \
    barvinok/evalue.h \
    barvinok/genfun.h \
    barvinok/options.h \
    barvinok/polylib.h \
    barvinok/sample.h \
    barvinok/basis_reduction.h
nodist_pkginclude_HEADERS = \
    barvinok/NTL.h \
    barvinok/set.h
LINK = $(CXXLINK)

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
if USE_ZSOLVE
    HILBERT_C = hilbert.c
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
    binomial.c \
    binomial.h \
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
    evalue_isl.c \
    genfun_constructor.cc \
    genfun_constructor.h \
    $(HILBERT_C) \
    hilbert.h \
    hull.c \
    hull.h \
    ilp.c \
    ilp.h \
    lattice_point.cc \
    lattice_point.h \
    lattice_width.c \
    lattice_width.h \
    laurent.cc \
    laurent.h \
    laurent_old.cc \
    laurent_old.h \
    normalization.c \
    normalization.h \
    options.c \
    param_polynomial.cc \
    param_polynomial.h \
    param_util.c \
    param_util.h \
    $(POLYSIGN_CDD) \
    $(POLYSIGN_GLPK) \
    polysign.c \
    polysign_pip.c \
    polysign_polylib.c \
    polysign.h \
    power.h \
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
    section_array.h \
    series.cc \
    $(TOPCOM) \
    summate.c \
    summate.h \
    vertex_cone.cc \
    vertex_cone.h \
    volume.c \
    volume.h \
    mat_util.cc \
    mat_util.h \
    matrix_read.cc \
    matrix_read.h \
    barvinok.cc \
    barvinok_e.cc \
    $(BR_CDD) \
    $(BR_GLPK)
EXTRA_libbarvinok_la_SOURCES = \
    basis_reduction_cdd.cc \
    basis_reduction_glpk.c \
    hilbert.c \
    polysign_cddf.cc \
    polysign_cdd.cc \
    polysign.c \
    initcdd.cc \
    initcdd.h \
    topcom.c \
    topcom.h
if BUNDLED_POLYLIB
    POLYLIB_LA = $(top_builddir)/polylib/libpolylibgmp.la
endif
if BUNDLED_PIPLIB
    PIPLIB_LA = $(top_builddir)/piplib/libpiplibMP.la
endif
if USE_ZSOLVE
    ZSOLVE_LA = zsolve/libzsolve.la
endif
libbarvinok_la_LIBADD = @LTLIBOBJS@ $(LIB_ISL_POLYLIB) $(LIB_ISL) \
			     $(POLYLIB_LA) @POLYLIB_LIBS@ \
			     $(PIPLIB_LA) @PIPLIB_LIBS@ \
			     $(ZSOLVE_LA)
libbarvinok_la_LDFLAGS = @BV_LDFLAGS@ -version-info @versioninfo@ \
				$(AM_LDFLAGS)
LDADD = libbarvinok.la @ISL_POLYLIB_LIBS@ @ISL_LIBS@ \
	@POLYLIB_LIBS@ @PIPLIB_LIBS@

test_SOURCES = \
	evalue_read.c \
	evalue_read.h \
	test.c
testlib_SOURCES = \
	evalue_read.c \
	evalue_read.h \
	testlib.cc
barvinok_count_SOURCES = barvinok_count.c
barvinok_ehrhart_SOURCES = \
	barvinok_ehrhart_options.c \
	barvinok_ehrhart_options.h \
	barvinok_ehrhart.cc \
	evalue_convert_options.c \
	evalue_convert.cc \
	evalue_convert.h \
	$(FDSTREAM)
barvinok_union_SOURCES = \
	barvinok_union_options.c \
	barvinok_union_options.h \
	barvinok_union.cc
if HAVE_OMEGA
BEEO_SOURCES = \
	omega_interface/Exit.cc \
	omega_interface/convert.h \
	omega_interface/convert.cc \
	omega_interface/count.h \
	omega_interface/count.cc
else
BEEO_SOURCES = 
endif
if USE_FDSTREAM
FDSTREAM = fdstream.cc fdstream.h
endif
barvinok_enumerate_SOURCES = \
	barvinok_enumerate_options.c \
	barvinok_enumerate_options.h \
	barvinok_enumerate.cc \
	evalue_convert_options.c \
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
	barvinok_enumerate_e_options.c \
	barvinok_enumerate_e_options.h \
	barvinok_enumerate_e.cc \
	evalue_convert_options.c \
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
if USE_PARKER
    MONA_LIBS = -ldfa -lbdd -lmem
    PARKER_LA = parker/libparker.la
    PARKERLDFLAGS = -L$(MONAPATH)/DFA -L$(MONAPATH)/BDD -L$(MONAPATH)/Mem
    PARKERCPPFLAGS = -I$(MONAPATH)/Mem -I$(MONAPATH)/DFA -I$(MONAPATH)/BDD
endif
barvinok_enumerate_e_CPPFLAGS = \
	$(PARKERCPPFLAGS) @OMEGA_CPPFLAGS@ $(AM_CPPFLAGS)
barvinok_enumerate_e_LDFLAGS = $(PARKERLDFLAGS) @OMEGA_LDFLAGS@ @LDFLAGS@
barvinok_enumerate_e_LDADD = \
	$(PARKER_LA) $(MONA_LIBS) @OMEGA_LIBS@ $(OMEGA_LA) libbarvinok.la \
	@ISL_POLYLIB_LIBS@ @ISL_LIBS@ @POLYLIB_LIBS@ @PIPLIB_LIBS@
barvinok_bound_SOURCES = \
	evalue_convert_options.c \
	evalue_convert.cc \
	evalue_convert.h \
	evalue_read.c \
	evalue_read.h \
	$(FDSTREAM) \
	bound_options.h \
	bound_options.c \
	bound.cc \
	verify.h \
	verify.c
EXTRA_barvinok_bound_SOURCES = fdstream.cc fdstream.h
barvinok_bound_LDADD = \
	libbarvinok.la @ISL_POLYLIB_LIBS@ @ISL_LIBS@ \
	@POLYLIB_LIBS@ @PIPLIB_LIBS@
barvinok_summate_SOURCES = \
	barvinok_summate_options.c \
	barvinok_summate_options.h \
	barvinok_summate.cc \
	evalue_convert_options.c \
	evalue_convert.cc \
	evalue_convert.h \
	evalue_read.c \
	evalue_read.h \
	$(FDSTREAM) \
	verify.h \
	verify.c
EXTRA_barvinok_summate_SOURCES = fdstream.cc fdstream.h
4coins_SOURCES = 4coins.cc
semigroup_holes_SOURCES = semigroup_holes.cc
lexmin_SOURCES = \
	lexmin_options.c \
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
test_bound_SOURCES = \
	test_bound.c \
	evalue_convert_options.c \
	evalue_convert.cc \
	evalue_convert.h \
	evalue_read.c \
	evalue_read.h \
	$(FDSTREAM) \
	verify.h \
	verify.c
EXTRA_test_bound_SOURCES = fdstream.cc fdstream.h
test_bound_LDADD = \
	libbarvinok.la @ISL_POLYLIB_LIBS@ @ISL_LIBS@ \
	@POLYLIB_LIBS@ @PIPLIB_LIBS@
evalue_convert_SOURCES = \
	evalue_convert_main.c \
	evalue_convert_options.c \
	evalue_convert.cc \
	evalue_convert.h \
	evalue_read.c \
	evalue_read.h \
	$(FDSTREAM)
iscc_SOURCES =	\
	isl_obj_list.h \
	isl_obj_list.c \
	iscc.c
iscc_LDADD = \
	libbarvinok.la @ISL_POLYLIB_LIBS@ @ISL_LIBS@ \
	@POLYLIB_LIBS@ @PIPLIB_LIBS@

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = \
    ChangeLog \
    $(TESTFILES) \
    latte2polylib.pl \
    NTL_5_3_2.patch \
    basis_reduction_templ.c \
    cdd94e-test \
    polysign_cdd_template.cc \
    barvinok/NTL.h.broken \
    barvinok/NTL.h.normal \
    barvinok/set.h.broken \
    barvinok/set.h.normal \
    doc/isl.tex \
    doc/Internal.tex \
    doc/applications.tex \
    doc/barvinok.tex \
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
	check-lexmin check-approx check-evalue check-euler check-lw
check-testlib: testlib$(EXEEXT)
	./testlib$(EXEEXT)
check-approx: test_approx$(EXEEXT)
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i | ./test_approx$(EXEEXT) -q -r 10 || exit; \
	    fi \
	done
check-evalue: barvinok_bound$(EXEEXT) test_bound$(EXEEXT) barvinok_summate$(EXEEXT)
	@for i in $(top_srcdir)/tests/evalue/*; do \
	    if test -f $$i; then \
		echo $$i; \
		if test -n "@bv_barvinok_bound@"; then \
		    ./barvinok_bound$(EXEEXT) -T -r30 < $$i || exit; \
		    echo $$i | ./test_bound$(EXEEXT) -q -r30 || exit; \
		fi; \
		./barvinok_summate$(EXEEXT) -T -r30 < $$i || exit; \
	    fi \
	done
check-euler: barvinok_summate$(EXEEXT)
	@for i in $(top_srcdir)/tests/euler/*; do \
	    if test -f $$i; then \
		for method in 'euler' 'laurent_old' 'laurent'; do \
		    opt="--summation=$$method"; \
		    echo $$i $$opt; \
		    ./barvinok_summate$(EXEEXT) -T $$opt < $$i || exit; \
		done \
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
		done; \
		opt="--summation=bernoulli"; \
		echo "        $$opt"; \
		./barvinok_enumerate$(EXEEXT) --verify $$opt < $$i || exit; \
	    fi \
	done
if HAVE_OMEGA
    SPEC_OMEGA = '--isl --omega'
endif
check-enumerate_e: barvinok_enumerate_e$(EXEEXT)
	@for dir in $(BEE_TESTDIRS); do \
	    for i in $(top_srcdir)/tests/$$dir/*; do \
		if test -f $$i; then \
		    for options in '' $(SPEC_OMEGA) '--isl'; do \
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
		done; \
		echo $$i --series; \
		./barvinok_enumerate_e$(EXEEXT) --verify --series < $$i || exit; \
	    fi \
	done
check-test: test$(EXEEXT)
	@failed=0; \
	for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		for options in '--index=10' '--primal --index=10'; do \
		    for spec in 'random' 'bf' 'df' 'todd'; do \
			opt="--specialization=$$spec $$options"; \
			echo -n $$i $$opt; \
			./test$(EXEEXT) $$opt < $$i; \
			if test "$$?" -ne "0"; then \
			    failed=`expr $$failed + 1`; \
			    echo " NOT ok"; \
			else \
			    echo " ok"; \
			fi; \
		    done \
		done \
	    fi \
	done; \
	for i in $(top_srcdir)/tests/sample/*; do \
	    echo -n $$i; \
	    ./test$(EXEEXT) < $$i; \
	    if test "$$?" -ne "0"; then \
		failed=`expr $$failed + 1`; \
		echo " NOT ok"; \
	    else \
		echo " ok"; \
	    fi; \
	done; \
	if test $$failed != 0; then \
		echo "$$failed tests failed"; \
		exit -1; \
	fi
check-lw: test$(EXEEXT)
	@failed=0; \
	for i in $(top_srcdir)/tests/lattice_width/*; do \
	    if test -f $$i; then \
		echo -n $$i; \
		./test$(EXEEXT) < $$i; \
		if test "$$?" -ne "0"; then \
		    failed=`expr $$failed + 1`; \
		    echo " NOT ok"; \
		else \
		    echo " ok"; \
		fi; \
	    fi \
	done; \
	if test $$failed != 0; then \
		echo "$$failed tests failed"; \
		exit -1; \
	fi
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
