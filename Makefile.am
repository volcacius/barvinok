lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count randomtest barvinok_enumerate \
		  verif_ehrhart barvinok_enumerate_e verif_ehrhart_e \
		  @bv_extra_programs@
EXTRA_PROGRAMS = piptest
pkginclude_HEADERS = barvinok.h util.h ev_operations.h
LINK = $(CXXLINK)

libbarvinok_la_SOURCES = \
    ev_operations.c \
    ev_operations.h \
    util.c \
    util.h \
    barvinok.cc \
    barvinok.h
EXTRA_libbarvinok_la_SOURCES = \
    piputil.h
libbarvinok_la_LIBADD = @LTLIBOBJS@
libbarvinok_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbarvinok.la

test_SOURCES = test.c
barvinok_count_SOURCES = barvinok_count.c

TESTFILES = $(shell find $(top_srcdir)/tests -type f ! -path '*CVS*')

EXTRA_DIST = ChangeLog $(TESTFILES)

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
