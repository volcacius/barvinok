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

bernstein_SOURCES = src/main.c
