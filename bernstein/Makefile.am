lib_LTLIBRARIES = libbernstein.la
bin_PROGRAMS = bernstein
pkginclude_HEADERS = include/bernstein/bernstein.h

AM_CPPFLAGS = -I$(srcdir)/include

libbernstein_la_SOURCES = \
    src/bernstein.cpp \
    src/piecewise_lst.cpp \
    src/piecewise_lst.h \
    src/polylib++.h \
    src/polynomial.cpp \
    src/polynomial.h
libbernstein_la_LIBADD = -lginac -lcln
libbernstein_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbernstein.la

bernstein_SOURCES = src/main.cpp

EXTRA_DIST = examples

check:
	@for i in $(top_srcdir)/examples/*.matrix; do \
	    if ./bernstein < $$i | diff -u $${i%.matrix}.out -; then \
		echo "$$i passed"; \
	    else \
		echo "$$i failed"; \
	    fi \
	done
