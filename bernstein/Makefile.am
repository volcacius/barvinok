lib_LTLIBRARIES = libbernstein.la
bin_PROGRAMS = bernstein
pkginclude_HEADERS = \
	include/bernstein/bernstein.h \
	include/bernstein/maximize.h \
	include/bernstein/piecewise_lst.h

AM_CPPFLAGS = -I$(srcdir)/include @GINACLIB_CPPFLAGS@

libbernstein_la_SOURCES = \
    src/bernstein.cpp \
    src/maximize.cpp \
    src/piecewise_lst.cpp \
    src/polylib++.h \
    src/polynomial.cpp \
    src/polynomial.h
libbernstein_la_LIBADD = @GINACLIB_LIBS@
libbernstein_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbernstein.la

bernstein_SOURCES = src/main.cpp

EXTRA_DIST = examples

ACLOCAL_AMFLAGS = -I m4

check:
	@for i in $(top_srcdir)/examples/*.matrix; do \
	    if ./bernstein < $$i | diff -u $${i%.matrix}.out -; then \
		echo "$$i passed"; \
	    else \
		echo "$$i failed"; \
	    fi \
	done
