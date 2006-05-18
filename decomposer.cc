#include <vector>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "decomposer.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif
using std::vector;

/*
 * Returns the largest absolute value in the vector
 */
static ZZ max(vec_ZZ& v)
{
    ZZ max = abs(v[0]);
    for (int i = 1; i < v.length(); ++i)
	if (abs(v[i]) > max)
	    max = abs(v[i]);
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
	matrix2zz(Rays, A, Rays->NbRows - 1, Rays->NbColumns - 1);
	det = determinant(A);
    }

    Vector* short_vector(vec_ZZ& lambda) {
	Matrix *M = Matrix_Copy(Rays);
	Matrix *inv = Matrix_Alloc(M->NbRows, M->NbColumns);
	int ok = Matrix_Inverse(M, inv);
	assert(ok);
	Matrix_Free(M);

	ZZ det2;
	mat_ZZ B;
	mat_ZZ U;
	matrix2zz(inv, B, inv->NbRows - 1, inv->NbColumns - 1);
	long r = LLL(det2, B, U);

	ZZ min = max(B[0]);
	int index = 0;
	for (int i = 1; i < B.NumRows(); ++i) {
	    ZZ tmp = max(B[i]);
	    if (tmp < min) {
		min = tmp;
		index = i;
	    }
	}

	Matrix_Free(inv);

	lambda = B[index];

	Vector *z = Vector_Alloc(U[index].length()+1);
	assert(z);
	zz2values(U[index], z->p);
	value_set_si(z->p[U[index].length()], 0);

	Polyhedron *C = poly();
	int i;
	for (i = 0; i < lambda.length(); ++i)
	    if (lambda[i] > 0)
		break;
	if (i == lambda.length()) {
	    Value tmp;
	    value_init(tmp);
	    value_set_si(tmp, -1);
	    Vector_Scale(z->p, z->p, tmp, z->Size-1);
	    value_clear(tmp);
	}
	return z;
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
	    Cone = Rays2Polyhedron(M, M->NbRows+1);
	    assert(Cone->NbConstraints == Cone->NbRays);
	    Matrix_Free(M);
	}
	return Cone;
    }

    ZZ det;
    Polyhedron *Cone;
    Matrix *Rays;
};

void decomposer::decompose(Polyhedron *C)
{
    vector<cone *> nonuni;
    cone * c = new cone(C);
    ZZ det = c->det;
    int s = sign(det);
    assert(det != 0);
    if (abs(det) > 1) {
	nonuni.push_back(c);
    } else {
	try {
	    handle(C, 1);
	    delete c;
	} catch (...) {
	    delete c;
	    throw;
	}
    }
    vec_ZZ lambda;
    while (!nonuni.empty()) {
	c = nonuni.back();
	nonuni.pop_back();
	Vector* v = c->short_vector(lambda);
	for (int i = 0; i < c->Rays->NbRows - 1; ++i) {
	    if (lambda[i] == 0)
		continue;
	    Matrix* M = Matrix_Copy(c->Rays);
	    Vector_Copy(v->p, M->p[i], v->Size);
	    cone * pc = new cone(M);
	    assert (pc->det != 0);
	    if (abs(pc->det) > 1) {
		assert(abs(pc->det) < abs(c->det));
		nonuni.push_back(pc);
	    } else {
		try {
		    handle(pc->poly(), sign(pc->det) * s);
		    delete pc;
		} catch (...) {
		    delete c;
		    delete pc;
		    while (!nonuni.empty()) {
			c = nonuni.back();
			nonuni.pop_back();
			delete c;
		    }
		    Matrix_Free(M);
		    Vector_Free(v);
		    throw;
		}
	    }
	    Matrix_Free(M);
	}
	Vector_Free(v);
	delete c;
    }
}

void polar_decomposer::decompose(Polyhedron *cone, unsigned MaxRays)
{
    Polyhedron_Polarize(cone);
    if (cone->NbRays - 1 != cone->Dimension) {
	Polyhedron *tmp = cone;
	cone = triangulate_cone(cone, MaxRays);
	Polyhedron_Free(tmp);
    }
    try {
	for (Polyhedron *Polar = cone; Polar; Polar = Polar->next)
	    decomposer::decompose(Polar);
	Domain_Free(cone);
    } catch (...) {
	Domain_Free(cone);
	throw;
    }
}

void polar_decomposer::handle(Polyhedron *P, int sign)
{
    Polyhedron_Polarize(P);
    handle_polar(P, sign);
}

void vertex_decomposer::decompose_at_vertex(Param_Vertices *V, int _i, 
					    unsigned MaxRays)
{
    Polyhedron *C = supporting_cone_p(P, V);
    vert = _i;
    this->V = V;

    pd->decompose(C, MaxRays);
}

/*
 * Barvinok's Decomposition of a simplicial cone
 *
 * Returns two lists of polyhedra
 */
void barvinok_decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg)
{
    Polyhedron *pos = *ppos, *neg = *pneg;
    vector<cone *> nonuni;
    cone * c = new cone(C);
    ZZ det = c->det;
    int s = sign(det);
    assert(det != 0);
    if (abs(det) > 1) {
	nonuni.push_back(c);
    } else {
	Polyhedron *p = Polyhedron_Copy(c->Cone);
	p->next = pos;
	pos = p;
	delete c;
    }
    vec_ZZ lambda;
    while (!nonuni.empty()) {
	c = nonuni.back();
	nonuni.pop_back();
	Vector* v = c->short_vector(lambda);
	for (int i = 0; i < c->Rays->NbRows - 1; ++i) {
	    if (lambda[i] == 0)
		continue;
	    Matrix* M = Matrix_Copy(c->Rays);
	    Vector_Copy(v->p, M->p[i], v->Size);
	    cone * pc = new cone(M);
	    assert (pc->det != 0);
	    if (abs(pc->det) > 1) {
		assert(abs(pc->det) < abs(c->det));
		nonuni.push_back(pc);
	    } else {
		Polyhedron *p = pc->poly();
		pc->Cone = 0;
		if (sign(pc->det) == s) {
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
    *ppos = pos;
    *pneg = neg;
}

