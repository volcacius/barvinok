#include <assert.h>
#include <iostream>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <barvinok.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

using namespace NTL;
using std::cout;
using std::endl;

#define ALLOC(p) (((long *) (p))[0])
#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

static void value2zz(Value v, ZZ& z)
{
    int sa = v[0]._mp_size;
    int abs_sa = sa < 0 ? -sa : sa;

    _ntl_gsetlength(&z.rep, abs_sa);
    mp_limb_t * adata = DATA(z.rep);
    for (int i = 0; i < abs_sa; ++i)
	adata[i] = v[0]._mp_d[i];
    SIZE(z.rep) = sa;
}

static void zz2value(ZZ& z, Value& v)
{
    int sa = SIZE(z.rep);
    int abs_sa = sa < 0 ? -sa : sa;

    mp_limb_t * adata = DATA(z.rep);
    mpz_realloc2(v, __GMP_BITS_PER_MP_LIMB * abs_sa);
    for (int i = 0; i < abs_sa; ++i)
	v[0]._mp_d[i] = adata[i];
    v[0]._mp_size = sa;
}

/*
 * We just ignore the last column and row
 * If the final element is not equal to one
 * then the result will actually be a multiple of the input
 */
static void matrix2zz(Matrix *M, mat_ZZ& m)
{
    m.SetDims(M->NbRows - 1, M->NbColumns - 1);

    for (int i = 0; i < M->NbRows-1; ++i) {
//	assert(value_one_p(M->p[i][M->NbColumns - 1]));
	for (int j = 0; j < M->NbColumns - 1; ++j) {
	    value2zz(M->p[i][j], m[i][j]);
	}
    }
}

static void rays(mat_ZZ& r, Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);
    r.SetDims(dim, dim);
    ZZ tmp;

    int i, c;
    for (i = 0, c = 0; i < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1])) {
	    for (int j = 0; j < dim; ++j) {
		value2zz(C->Ray[i][j+1], tmp);
		r[j][c] = tmp;
	    }
	    ++c;
	}
}

static Matrix * rays(Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);

    Matrix *M = Matrix_Alloc(dim+1, dim+1);
    assert(M);

    int i, c;
    for (i = 0, c = 0; i < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1])) {
	    Vector_Copy(C->Ray[i] + 1, M->p[c], dim);
	    value_set_si(M->p[c++][dim], 1);
	}
    value_set_si(M->p[dim][dim], 1);

    return M;
}

static ZZ max(vec_ZZ& v)
{
    ZZ max = v[0];
    for (int i = 1; i < v.length(); ++i)
	if (v[i] > max)
	    max = v[i];
    return max;
}

class cone {
public:
    cone(Polyhedron *C) {
	Cone = Polyhedron_Copy(C);
	Rays = rays(C);
	Matrix_Print(stdout, P_VALUE_FMT, Rays);
	mat_ZZ A;
	matrix2zz(Rays, A);
	cout << A << endl;
	det = determinant(A);
	cout << "det: " << det << endl;
	Value v;
	value_init(v);
	zz2value(det, v);
	value_print(stdout, P_VALUE_FMT, v);
	puts("");
	value_clear(v);
    }

    void short_vector() {
	Matrix *M = Matrix_Copy(Rays);
	Matrix *inv = Matrix_Alloc(M->NbRows, M->NbColumns);
	int ok = Matrix_Inverse(M, inv);
	assert(ok);
	Matrix_Print(stdout, P_VALUE_FMT, inv);
	Matrix_Free(M);

	ZZ det2;
	mat_ZZ B;
	mat_ZZ U;
	matrix2zz(inv, B);
	long r = LLL(det2, B, U);

	cout << det2 << B << U << endl;

	ZZ min = max(U[0]);
	int index = 0;
	for (int i = 1; i < U.NumRows(); ++i) {
	    ZZ tmp = max(U[1]);
	    if (tmp < min) {
		min = tmp;
		index = i;
	    }
	}
	cout << index << ": " << min << endl;

	Matrix_Free(inv);
    }

    ~cone() {
	Polyhedron_Free(Cone);
	Matrix_Free(Rays);
    }

    ZZ det;
    Polyhedron *Cone;
    Matrix *Rays;
};

/*
 * Barvinok's Decomposition of a simplicial cone
 *
 * Returns a list of polyhedra
 */
Polyhedron *decompose(Polyhedron *C)
{
    mat_ZZ r;
    cone * c = new cone(C);
    c->short_vector();
    delete c;
}
