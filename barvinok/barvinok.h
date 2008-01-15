#ifndef BARVINOK_H
#define BARVINOK_H

#include <barvinok/evalue.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <barvinok/options.h>

void barvinok_decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg);
void barvinok_count(Polyhedron *P, Value* result, unsigned NbMaxCons);
void barvinok_count_with_options(Polyhedron *P, Value* result,
				 struct barvinok_options *options);
Enumeration* barvinok_enumerate(Polyhedron *P, Polyhedron* C, unsigned MaxRays);
evalue* barvinok_enumerate_ev(Polyhedron *P, Polyhedron* C, unsigned MaxRays);
evalue* barvinok_enumerate_with_options(Polyhedron *P, Polyhedron* C,
					struct barvinok_options *options);
evalue *Param_Polyhedron_Enumerate(Param_Polyhedron *PP, Polyhedron *P,
				   Polyhedron *C,
				   struct barvinok_options *options);
evalue* barvinok_enumerate_e(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays);
evalue* barvinok_enumerate_e_with_options(Polyhedron *P, 
		  unsigned exist, unsigned nparam, struct barvinok_options *options);
evalue *barvinok_enumerate_pip(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays);
evalue *barvinok_enumerate_pip_with_options(Polyhedron *P,
		  unsigned exist, unsigned nparam, struct barvinok_options *options);
evalue *barvinok_enumerate_scarf(Polyhedron *P,
			  unsigned exist, unsigned nparam,
			  struct barvinok_options *options);
evalue* barvinok_enumerate_union(Polyhedron *D, Polyhedron* C, unsigned MaxRays);

evalue *barvinok_summate(evalue *e, int nvar, struct barvinok_options *options);

/* private function */
evalue* bv_ceil3(Value *coef, int len, Value d, Polyhedron *P);

#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)

#include <barvinok/genfun.h>

void zz2value(const ZZ& z, Value& v);
gen_fun * barvinok_series(Polyhedron *P, Polyhedron* C, unsigned MaxRays);
gen_fun * barvinok_series_with_options(Polyhedron *P, Polyhedron* C,
				       barvinok_options *options);
gen_fun *barvinok_enumerate_series(Polyhedron *P, unsigned nparam,
				    barvinok_options *options);
gen_fun *barvinok_enumerate_scarf_series(Polyhedron *P,
			  unsigned exist, unsigned nparam,
			  barvinok_options *options);
gen_fun *barvinok_enumerate_e_series(Polyhedron *P,
		  unsigned exist, unsigned nparam, barvinok_options *options);
gen_fun* barvinok_enumerate_union_series(Polyhedron *D, Polyhedron* C, 
					 unsigned MaxRays);
gen_fun* barvinok_enumerate_union_series_with_options(Polyhedron *D, Polyhedron* C, 
						      barvinok_options *options);

#endif

#endif
