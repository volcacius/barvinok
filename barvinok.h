#ifndef BARVINOK_H
#define BARVINOK_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

void barvinok_decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg);
void barvinok_count(Polyhedron *P, Value* result, unsigned NbMaxCons);
Enumeration* barvinok_enumerate(Polyhedron *P, Polyhedron* C, unsigned MaxRays);
evalue* barvinok_enumerate_e(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays);

#ifdef USE_MODULO
evalue* ceil3(Value *coef, int len, Value d);
#endif

#if defined(__cplusplus)
}
#endif

#endif
