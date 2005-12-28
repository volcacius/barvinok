#ifndef BARVINOK_H
#define BARVINOK_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <barvinok/evalue.h>

void barvinok_decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg);
void barvinok_count(Polyhedron *P, Value* result, unsigned NbMaxCons);
Enumeration* barvinok_enumerate(Polyhedron *P, Polyhedron* C, unsigned MaxRays);
evalue* barvinok_enumerate_ev(Polyhedron *P, Polyhedron* C, unsigned MaxRays);
evalue* barvinok_enumerate_e(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays);
evalue *barvinok_enumerate_pip(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays);

/* private function */
evalue* bv_ceil3(Value *coef, int len, Value d, Polyhedron *P);

#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)

#include <barvinok/genfun.h>

void zz2value(ZZ& z, Value& v);
gen_fun * barvinok_series(Polyhedron *P, Polyhedron* C, unsigned MaxRays);

#endif

#endif
