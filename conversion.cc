#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include "conversion.h"

#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

void value2zz(Value v, ZZ& z)
{
    int sa = v[0]._mp_size;
    int abs_sa = sa < 0 ? -sa : sa;

    _ntl_gsetlength(&z.rep, abs_sa);
    mp_limb_t * adata = DATA(z.rep);
    for (int i = 0; i < abs_sa; ++i)
	adata[i] = v[0]._mp_d[i];
    SIZE(z.rep) = sa;
}

void zz2value(ZZ& z, Value& v)
{
    if (!z.rep) {
	value_set_si(v, 0);
	return;
    }

    int sa = SIZE(z.rep);
    int abs_sa = sa < 0 ? -sa : sa;

    mp_limb_t * adata = DATA(z.rep);
    _mpz_realloc(v, abs_sa);
    for (int i = 0; i < abs_sa; ++i)
	v[0]._mp_d[i] = adata[i];
    v[0]._mp_size = sa;
}

void values2zz(Value *p, vec_ZZ& v, int len)
{
    v.SetLength(len);

    for (int i = 0; i < len; ++i) {
	value2zz(p[i], v[i]);
    }
}

/*
 */
void zz2values(vec_ZZ& v, Value *p)
{
    for (int i = 0; i < v.length(); ++i)
	zz2value(v[i], p[i]);
}

/*
 * We just ignore the last column and row
 * If the final element is not equal to one
 * then the result will actually be a multiple of the input
 */
void matrix2zz(Matrix *M, mat_ZZ& m, unsigned nr, unsigned nc)
{
    m.SetDims(nr, nc);

    for (int i = 0; i < nr; ++i) {
//	assert(value_one_p(M->p[i][M->NbColumns - 1]));
	for (int j = 0; j < nc; ++j) {
	    value2zz(M->p[i][j], m[i][j]);
	}
    }
}

Matrix *rays(Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);

    Matrix *M = Matrix_Alloc(dim+1, dim+1);
    assert(M);

    int i, c;
    for (i = 0, c = 0; i <= dim && c < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1])) {
	    Vector_Copy(C->Ray[i] + 1, M->p[c], dim);
	    value_set_si(M->p[c++][dim], 0);
	}
    assert(c == dim);
    value_set_si(M->p[dim][dim], 1);

    return M;
}

Matrix * rays2(Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);

    Matrix *M = Matrix_Alloc(dim, dim);
    assert(M);

    int i, c;
    for (i = 0, c = 0; i <= dim && c < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1]))
	    Vector_Copy(C->Ray[i] + 1, M->p[c++], dim);
    assert(c == dim);

    return M;
}
