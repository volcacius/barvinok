# ------------------------------------------------
#
#	Bernstein Expansion
#
# ------------------------------------------------


# Debug flag
#DEBUG=-DDEBUG=1

# Sets for the polylib
POLYLIB = /usr/local
POLYLIB_SOFT = /usr/local/soft/polylib
POLYLIB_FLAGS = -DSIZEOF_INT=4 -DSIZEOF_LONG_INT=4 -DSIZEOF_LONG_LONG_INT=8 -DGNUMP

# Sets for gmp
GMP_SOFT = /usr/local/soft/gmp
GMP = /usr

# Include paths
INCLUDES = -I$(POLYLIB)/include -I$(POLYLIB_SOFT)/include -I$(GMP)/include -I$(GMP_SOFT)/include

# Libs directory
LIBS= -L$(POLYLIB)/lib -L$(GMP)/lib -L/usr/local/lib -L$(POLYLIB_SOFT)/lib -L$(GMP_SOFT)/lib

# Libs to include when linking
LIBS_TO_ADD = -lginac -lcln -lpolylibgmp -lgmp

# Flags for the compiler
CFLAGS = $(POLYLIB_FLAGS) $(DEBUG)

#flags for the linker
LFLAGS = 

# program name
PROGNAME=bernstein

# Compilers
CC=gcc
CPPC=g++

# Warnings
WARNINGS=-Wall

# Object files
OBJS=bernstein.o bernstein++.o bernstein-expansion.o polynomial.o

# quizas seria bueno poner esto para la cnl (en mi caso) `-Wl,--rpath -Wl,LIBDIR'


all: $(PROGNAME)

$(PROGNAME): $(OBJS)
	$(CPPC) $(OBJS) -o $@ $(CFLAGS) $(LIBS_TO_ADD) $(LIBS) $(LFLAGS)
clean:
	rm -f core $(PROGNAME) $(OBJS) *~

strip:
	strip $(PROGNAME)

.cpp.o:
	$(CPPC) $(CFLAGS) $(WARNINGS) $(INCLUDES) -c $<

.c.o:
	$(CC) $(CFLAGS) $(WARNINGS) $(INCLUDES) -c $<
