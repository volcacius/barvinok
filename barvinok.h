#ifndef BARVINOK_H
#define BARVINOK_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

void decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg);
void count(Polyhedron *P);

#if defined(__cplusplus)
}
#endif

#endif
