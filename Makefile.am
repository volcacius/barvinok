lib_LTLIBRARIES = libbarvinok.la
noinst_PROGRAMS = test barvinok_count

libbarvinok_la_SOURCES = \
    util.c \
    util.h \
    barvinok.cc \
    barvinok.h
libbarvinok_la_LIBADD = @LTLIBOBJS@ -lstdc++
libbarvinok_la_LDFLAGS = -version-info @versioninfo@
LDADD = libbarvinok.la

test_SOURCES = test.c
barvinok_count_SOURCES = barvinok_count.c
