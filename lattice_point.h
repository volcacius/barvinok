#include <gmp.h>
#include <NTL/vec_ZZ.h>
#include <barvinok/polylib.h>

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

evalue *multi_monom(vec_ZZ& p);
int normal_mod(Value *coef, int len, Value *m);
void lattice_point(Value* values, const mat_ZZ& rays, vec_ZZ& vertex, int *closed);
evalue* lattice_point(const mat_ZZ& rays, vec_ZZ& lambda, Matrix *W,
		      Value lcm, Polyhedron *PD, barvinok_options *options);
void lattice_point(Param_Vertices *V, const mat_ZZ& rays, vec_ZZ& num, 
		   evalue **E_vertex, barvinok_options *options);
