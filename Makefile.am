lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count randomtest barvinok_enumerate

libbarvinok_la_SOURCES = \
    debug.cpp \
    debug.h \
    ehrhartpolynom.cpp \
    ehrhartpolynom.h \
    util.c \
    util.h \
    barvinok.cc \
    barvinok.h
libbarvinok_la_LIBADD = @LTLIBOBJS@ -lstdc++
libbarvinok_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbarvinok.la

test_SOURCES = test.c
barvinok_count_SOURCES = barvinok_count.c

TESTFILES = \
	$(filter-out $(top_srcdir)/tests/CVS,$(wildcard $(top_srcdir)/tests/*))

EXTRA_DIST = ChangeLog $(TESTFILES)

check:
	@for i in $(top_srcdir)/tests/*; do \
	    if test -f $$i; then \
		echo -n $$i; \
		./test < $$i || echo -n " NOT"; \
		echo " ok"; \
	    fi \
	done
