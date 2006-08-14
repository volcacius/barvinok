SUBDIRS = . omega doc

lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count randomtest barvinok_enumerate \
		  barvinok_ehrhart \
		  verif_ehrhart barvinok_enumerate_e \
		  barvinok_series remove_redundant_equalities \
		  barvinok_union \
		  @bv_extra_programs@
EXTRA_PROGRAMS = piptest verify_lexsmaller polyhedron_sample 4coins lexmin
pkginclude_HEADERS = \
    barvinok/NTL_QQ.h \
    barvinok/barvinok.h \
    barvinok/util.h \
    barvinok/evalue.h \
    barvinok/genfun.h
LINK = $(CXXLINK)

libbarvinok_la_SOURCES = \
    barvinok/evalue.h \
    barvinok/genfun.h \
    barvinok/util.h \
    barvinok/barvinok.h \
    NTL_QQ.cc \
    evalue.c \
    genfun.cc \
    util.c \
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
    reduce_domain.c \
    reduce_domain.h \
    reducer.cc \
    reducer.h \
    mat_util.cc \
    mat_util.h \
    barvinok.cc
EXTRA_libbarvinok_la_SOURCES = \
    piputil.h
libbarvinok_la_LIBADD = @LTLIBOBJS@
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
	barvinok_enumerate_e.cc \
	scarf.cc \
	$(BEEO_SOURCES)
barvinok_enumerate_e_CPPFLAGS = @OMEGA_CPPFLAGS@ @CPPFLAGS@
barvinok_enumerate_e_LDFLAGS = @OMEGA_LDFLAGS@ @LDFLAGS@
barvinok_enumerate_e_LDADD = @OMEGA_LIBS@ libbarvinok.la
verif_ehrhart_SOURCES = verif_ehrhart.c verify_main.cc verif_ehrhart.h
4coins_SOURCES = 4coins.cc scarf.cc
lexmin_SOURCES = lexmin.cc combine.c combine.h fdstream.cc fdstream.h

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = \
    README.Solaris \
    ChangeLog \
    $(TESTFILES) \
    latte2polylib.pl \
    NTL_5_3_2.patch \
    doc/Internal.tex \
    doc/Usage.tex \
    doc/barvinok.bib \
    doc/barvinok.gdf \
    doc/barvinok.tex \
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
BEE_TESTDIRS = ehrhart_e ehrhart_e/piplib
else
BEE_TESTDIRS = ehrhart_e
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
	@for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		echo -n $$i; \
		./test < $$i || echo -n " NOT"; \
		echo " ok"; \
	    fi \
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
