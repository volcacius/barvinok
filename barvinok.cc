#include <assert.h>
#include <iostream>
#include <vector>
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
using std::vector;

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
    if (!z.rep) {
	value_set_si(v, 0);
	return;
    }

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

/*
 * We add a 0 at the end, because we need it afterwards
 */
static Vector * zz2vector(vec_ZZ& v)
{
    Vector *vec = Vector_Alloc(v.length()+1);
    assert(vec);
    for (int i = 0; i < v.length(); ++i)
	zz2value(v[i], vec->p[i]);

    value_set_si(vec->p[v.length()], 0);

    return vec;
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
	    value_set_si(M->p[c++][dim], 0);
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
    cone(Matrix *M) {
	Cone = 0;
	Rays = Matrix_Copy(M);
	set_det();
    }
    cone(Polyhedron *C) {
	Cone = Polyhedron_Copy(C);
	Rays = rays(C);
	set_det();
    }
    void set_det() {
	mat_ZZ A;
	matrix2zz(Rays, A);
	det = determinant(A);
	Value v;
	value_init(v);
	zz2value(det, v);
	value_clear(v);
    }

    Vector* short_vector() {
	Matrix *M = Matrix_Copy(Rays);
	Matrix *inv = Matrix_Alloc(M->NbRows, M->NbColumns);
	int ok = Matrix_Inverse(M, inv);
	assert(ok);
	Matrix_Free(M);

	ZZ det2;
	mat_ZZ B;
	mat_ZZ U;
	matrix2zz(inv, B);
	long r = LLL(det2, B, U);

	ZZ min = max(U[0]);
	int index = 0;
	for (int i = 1; i < U.NumRows(); ++i) {
	    ZZ tmp = max(U[1]);
	    if (tmp < min) {
		min = tmp;
		index = i;
	    }
	}

	Matrix_Free(inv);
	return zz2vector(U[index]);
    }

    ~cone() {
	Polyhedron_Free(Cone);
	Matrix_Free(Rays);
    }

    Polyhedron *poly() {
	if (!Cone) {
	    Matrix *M = Matrix_Alloc(Rays->NbRows+1, Rays->NbColumns+1);
	    for (int i = 0; i < Rays->NbRows; ++i) {
		Vector_Copy(Rays->p[i], M->p[i]+1, Rays->NbColumns);
		value_set_si(M->p[i][0], 1);
	    }
	    Vector_Set(M->p[Rays->NbRows]+1, 0, Rays->NbColumns-1);
	    value_set_si(M->p[Rays->NbRows][0], 1);
	    value_set_si(M->p[Rays->NbRows][Rays->NbColumns], 1);
	    Cone = Rays2Polyhedron(M, 600);
	    Matrix_Free(M);
	}
	return Cone;
    }

    ZZ det;
    Polyhedron *Cone;
    Matrix *Rays;
};

/*
 * Barvinok's Decomposition of a simplicial cone
 *
 * Returns two lists of polyhedra
 */
void decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg)
{
    Polyhedron *pos = 0, *neg = 0;
    vector<cone *> nonuni;
    cone * c = new cone(C);
    ZZ det = c->det;
    if (abs(det) > 1) {
	nonuni.push_back(c);
    } else {
	if (det > 0) 
	    pos = Polyhedron_Copy(c->Cone);
	else
	    neg = Polyhedron_Copy(c->Cone);
	delete c;
    }
    while (!nonuni.empty()) {
	c = nonuni.back();
	nonuni.pop_back();
	Vector* v = c->short_vector();
	for (int i = 0; i < c->Rays->NbRows - 1; ++i) {
	    Matrix* M = Matrix_Copy(c->Rays);
	    Vector_Copy(v->p, M->p[i], v->Size);
	    cone * pc = new cone(M);
	    if (abs(pc->det) > 1)
		nonuni.push_back(pc);
	    else {
		Polyhedron *p = Polyhedron_Copy(pc->poly());
		if (pc->det > 0) {
		    p->next = pos;
		    pos = p;
		} else {
		    p->next = neg;
		    neg = p;
		}
		delete pc;
	    }
	    Matrix_Free(M);
	}
	Vector_Free(v);
	delete c;
    }
    if (det > 0) {
	*ppos = pos;
	*pneg = neg;
    } else {
	*ppos = neg;
	*pneg = pos;
    }
}
