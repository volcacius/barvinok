lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count randomtest barvinok_enumerate \
		  verif_ehrhart barvinok_enumerate_e verif_ehrhart_e \
		  barvinok_series remove_redundant_equalities \
		  @bv_extra_programs@
EXTRA_PROGRAMS = piptest
pkginclude_HEADERS = barvinok.h util.h ev_operations.h
LINK = $(CXXLINK)

libbarvinok_la_SOURCES = \
    ev_operations.c \
    ev_operations.h \
    genfun.cc \
    genfun.h \
    util.c \
    util.h \
    barvinok.cc \
    barvinok.h \
    barvinok2.h
EXTRA_libbarvinok_la_SOURCES = \
    piputil.h
libbarvinok_la_LIBADD = @LTLIBOBJS@
libbarvinok_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbarvinok.la

test_SOURCES = test.c
barvinok_count_SOURCES = barvinok_count.c
barvinok_series_SOURCES = barvinok_series.cc

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = ChangeLog $(TESTFILES) cdd2polylib.pl NTL_5_3_2.patch

check:
	@for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		echo -n $$i; \
		./test < $$i || echo -n " NOT"; \
		echo " ok"; \
	    fi \
	done
	@for i in $(top_srcdir)/tests/ehrhart/*; do \
	    if test -f $$i; then \
		echo $$i; \
		./verif_ehrhart < $$i; \
	    fi \
	done
