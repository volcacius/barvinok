#include <gmp.h>
#include <NTL/vec_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

evalue *multi_monom(vec_ZZ& p);
int normal_mod(Value *coef, int len, Value *m);
void lattice_point(Value* values, Polyhedron *i, vec_ZZ& vertex);
evalue* lattice_point(
    Polyhedron *i, vec_ZZ& lambda, Matrix *W, Value lcm, Polyhedron *PD);
void lattice_point(Param_Vertices *V, Polyhedron *C, vec_ZZ& num, 
		   evalue **E_vertex);
