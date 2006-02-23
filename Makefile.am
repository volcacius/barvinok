lib_LTLIBRARIES = libbernstein.la
bin_PROGRAMS = bernstein

libbernstein_la_SOURCES = \
    src/bernstein++.cpp \
    src/bernstein++.h \
    src/bernstein-expansion.cpp \
    src/bernstein-expansion.h \
    src/bernstein.c \
    src/bernstein.h \
    src/polynomial.cpp \
    src/polynomial.h
libbernstein_la_LIBADD = -lginac -lcln
libbernstein_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbernstein.la

bernstein_SOURCES = src/main.cpp

check:
	@for i in $(top_srcdir)/examples/*.matrix; do \
	    if ./bernstein < $$i | diff -u - $${i%.matrix}.out; then \
		echo "$$i passed"; \
	    else \
		echo "$$i failed"; \
	    fi \
	done
