lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count randomtest barvinok_enumerate \
		  verif_ehrhart barvinok_enumerate_e verif_ehrhart_e \
		  barvinok_series remove_redundant_equalities \
		  verify_ranking \
		  @bv_extra_programs@
EXTRA_PROGRAMS = piptest
pkginclude_HEADERS = \
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
    evalue.c \
    genfun.cc \
    util.c \
    barvinok.cc
EXTRA_libbarvinok_la_SOURCES = \
    piputil.h
libbarvinok_la_LIBADD = @LTLIBOBJS@
libbarvinok_la_LDFLAGS = @BV_LDFLAGS@ -version-info @versioninfo@
LDADD = libbarvinok.la

test_SOURCES = test.c
barvinok_count_SOURCES = barvinok_count.c
barvinok_series_SOURCES = barvinok_series.cc
verif_ehrhart_SOURCES = verif_ehrhart.c verify_main.cc verif_ehrhart.h

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = \
    ChangeLog \
    $(TESTFILES) \
    cdd2polylib.pl \
    NTL_5_3_2.patch \
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
	if test -f $(top_srcdir)/.git/HEAD; then \
	    cp $(top_srcdir)/.git/HEAD $(distdir)/GIT_HEAD_ID; \
	else \
	    echo $(GIT_HEAD_ID) > $(distdir)/GIT_HEAD_ID; \
	fi

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
	@for i in $(top_srcdir)/tests/ehrhart_e/*; do \
	    if test -f $$i; then \
		echo $$i; \
		./verif_ehrhart_e < $$i; \
	    fi \
	done
	@for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		echo -n $$i; \
		./test < $$i || echo -n " NOT"; \
		echo " ok"; \
	    fi \
	done

version.h: @GIT_HEAD@
	echo '#define GIT_HEAD_ID "'@GIT_HEAD_VERSION@'"' > $@
