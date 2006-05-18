#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

void value2zz(Value v, ZZ& z);
void zz2value(ZZ& z, Value& v);
void zz2values(vec_ZZ& v, Value *p);
void matrix2zz(Matrix *M, mat_ZZ& m, unsigned nr, unsigned nc);
Matrix *rays(Polyhedron *C);
