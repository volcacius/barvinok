#include <vector>
#include <assert.h>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "decomposer.h"
#include "reduce_domain.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif
using std::vector;
using std::cerr;
using std::endl;

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

/* Remove common divisor of elements of cols of B */
static void normalize_cols(mat_ZZ& B)
{
    ZZ gcd;
    for (int i = 0; i < B.NumCols(); ++i) {
	gcd = B[0][i];
	for (int j = 1 ; gcd != 1 && j < B.NumRows(); ++j)
	    GCD(gcd, gcd, B[j][i]);
	if (gcd != 1)
	    for (int j = 0; j < B.NumRows(); ++j)
		B[j][i] /= gcd;
    }
}

class cone {
public:
    cone(const mat_ZZ& r, int row, const vec_ZZ& w) {
	rays = r;
	rays[row] = w;
	set_det();
	set_closed(NULL);
    }
    cone(const signed_cone& sc) {
	rays = sc.rays;
	set_det();
	set_closed(sc.closed);
    }
    void set_det() {
	det = determinant(rays);
    }
    void set_closed(int *cl) {
	closed = NULL;
	if (cl) {
	    closed = new int[rays.NumRows()];
	    for (int i = 0; i < rays.NumRows(); ++i)
		closed[i] = cl[i];
	}
    }
    bool needs_split(barvinok_options *options) {
	index = abs(det);
	if (IsOne(index))
	    return false;
	if (options->primal && index <= options->max_index)
	    return false;

	Matrix *M = zz2matrix(rays);
	Matrix *inv = Matrix_Alloc(M->NbRows, M->NbColumns);
	int ok = Matrix_Inverse(M, inv);
	assert(ok);
	Matrix_Free(M);

	matrix2zz(inv, B, inv->NbRows, inv->NbColumns);
	Matrix_Free(inv);

	if (!options->primal && options->max_index > 1) {
	    mat_ZZ B2 = B;
	    normalize_cols(B2);
	    index = abs(determinant(B2));
	    if (index <= options->max_index)
		return false;
	}

	return true;
    }

    void short_vector(vec_ZZ& v, vec_ZZ& lambda, barvinok_options *options) {
	ZZ det2;
	mat_ZZ U;
	long r = LLL(det2, B, U, options->LLL_a, options->LLL_b);

	ZZ min = max(B[0]);
	int index = 0;
	for (int i = 1; i < B.NumRows(); ++i) {
	    ZZ tmp = max(B[i]);
	    if (tmp < min) {
		min = tmp;
		index = i;
	    }
	}

	lambda = B[index];

	v = U[index];

	int i;
	for (i = 0; i < lambda.length(); ++i)
	    if (lambda[i] > 0)
		break;
	if (i == lambda.length()) {
	    v = -v;
	    lambda = -lambda;
	}
    }

    ~cone() {
	if (closed)
	    delete [] closed;
    }

    ZZ det;
    ZZ index;
    mat_ZZ rays;
    mat_ZZ B;
    int *closed;
};

static void decompose(const signed_cone& sc, signed_cone_consumer& scc,
			     barvinok_options *options)
{
    vector<cone *> nonuni;
    cone * c = new cone(sc);
    ZZ det = c->det;
    int s = sign(det);
    assert(det != 0);
    if (c->needs_split(options)) {
	nonuni.push_back(c);
    } else {
	try {
	    options->stats->base_cones++;
	    scc.handle(signed_cone(sc.C, sc.rays, sc.sign, to_ulong(c->index),
				   sc.closed), options);
	    delete c;
	} catch (...) {
	    delete c;
	    throw;
	}
	return;
    }
    vec_ZZ lambda;
    vec_ZZ v;;
    int closed[c->rays.NumRows()];
    while (!nonuni.empty()) {
	c = nonuni.back();
	nonuni.pop_back();
	c->short_vector(v, lambda, options);
	for (int i = 0; i < c->rays.NumRows(); ++i) {
	    if (lambda[i] == 0)
		continue;
	    cone *pc = new cone(c->rays, i, v);
	    if (c->closed) {
		for (int j = 0; j < c->rays.NumRows(); ++j) {
		    if (lambda[j] == 0)
			closed[j] = c->closed[j];
		    else if (j == i) {
			if (lambda[i] > 0)
			    closed[j] = c->closed[j];
			else
			    closed[j] = !c->closed[j];
		    } else if (sign(lambda[i]) == sign(lambda[j])) {
			if (c->closed[i] == c->closed[j])
			    closed[j] = i < j;
			else
			    closed[j] = c->closed[j];
		    } else
			closed[j] = c->closed[i] && c->closed[j];
		}
		pc->set_closed(closed);
	    }
	    assert (pc->det != 0);
	    if (pc->needs_split(options)) {
		assert(abs(pc->det) < abs(c->det));
		nonuni.push_back(pc);
	    } else {
		try {
		    options->stats->base_cones++;
		    scc.handle(signed_cone(pc->rays, sign(pc->det) * s,
					   to_ulong(pc->index),
					   pc->closed), options);
		    delete pc;
		} catch (...) {
		    delete c;
		    delete pc;
		    while (!nonuni.empty()) {
			c = nonuni.back();
			nonuni.pop_back();
			delete c;
		    }
		    throw;
		}
	    }
	}
	delete c;
    }
}

struct polar_signed_cone_consumer : public signed_cone_consumer {
    signed_cone_consumer& scc;
    mat_ZZ r;
    polar_signed_cone_consumer(signed_cone_consumer& scc) : scc(scc) {}
    virtual void handle(const signed_cone& sc, barvinok_options *options) {
	Polyhedron *C = sc.C;
	if (!sc.C) {
	    Matrix *M = Matrix_Alloc(sc.rays.NumRows()+1, sc.rays.NumCols()+2);
	    for (int i = 0; i < sc.rays.NumRows(); ++i) {
		zz2values(sc.rays[i], M->p[i]+1);
		value_set_si(M->p[i][0], 1);
	    }
	    value_set_si(M->p[sc.rays.NumRows()][0], 1);
	    value_set_si(M->p[sc.rays.NumRows()][1+sc.rays.NumCols()], 1);
	    C = Rays2Polyhedron(M, M->NbRows+1);
	    assert(C->NbConstraints == C->NbRays);
	    Matrix_Free(M);
	}
	Polyhedron_Polarize(C);
	rays(C, r);
	scc.handle(signed_cone(C, r, sc.sign, sc.det), options);
	if (!sc.C)
	    Polyhedron_Free(C);
    }
};

/* Remove common divisor of elements of rows of B */
static void normalize_rows(mat_ZZ& B)
{
    ZZ gcd;
    for (int i = 0; i < B.NumRows(); ++i) {
	gcd = B[i][0];
	for (int j = 1 ; gcd != 1 && j < B.NumCols(); ++j)
	    GCD(gcd, gcd, B[i][j]);
	if (gcd != 1)
	    for (int j = 0; j < B.NumCols(); ++j)
		B[i][j] /= gcd;
    }
}

static void polar_decompose(Polyhedron *cone, signed_cone_consumer& scc,
			    barvinok_options *options)
{
    POL_ENSURE_VERTICES(cone);
    Polyhedron_Polarize(cone);
    if (cone->NbRays - 1 != cone->Dimension) {
	Polyhedron *tmp = cone;
	cone = triangulate_cone_with_options(cone, options);
	Polyhedron_Free(tmp);
    }
    polar_signed_cone_consumer pssc(scc);
    mat_ZZ r;
    try {
	for (Polyhedron *Polar = cone; Polar; Polar = Polar->next) {
	    rays(Polar, r);
	    normalize_rows(r);
	    decompose(signed_cone(Polar, r, 1), pssc, options);
	}
	Domain_Free(cone);
    } catch (...) {
	Domain_Free(cone);
	throw;
    }
}

static void primal_decompose(Polyhedron *cone, signed_cone_consumer& scc,
			     barvinok_options *options)
{
    POL_ENSURE_VERTICES(cone);
    Polyhedron *parts;
    if (cone->NbRays - 1 == cone->Dimension)
	parts = cone;
    else
	parts = triangulate_cone_with_options(cone, options);
    int closed[cone->Dimension];
    Vector *average = NULL;
    Value tmp;
    if (parts != cone) {
	value_init(tmp);
	average = inner_point(cone);
    }
    mat_ZZ ray;
    try {
	for (Polyhedron *simple = parts; simple; simple = simple->next) {
	    for (int i = 0, r = 0; r < simple->NbRays; ++r) {
		if (value_notzero_p(simple->Ray[r][1+simple->Dimension]))
		    continue;
		if (simple == cone) {
		    closed[i] = 1;
		} else {
		    int f;
		    for (f = 0; f < simple->NbConstraints; ++f) {
			Inner_Product(simple->Ray[r]+1, simple->Constraint[f]+1,
				      simple->Dimension, &tmp);
			if (value_notzero_p(tmp))
			    break;
		    }
		    assert(f < simple->NbConstraints);
		    closed[i] = is_internal(average, simple->Constraint[f]);
		}
		++i;
	    }
	    rays(simple, ray);
	    decompose(signed_cone(simple, ray, 1, 0, closed), scc, options);
	}
	Domain_Free(parts);
	if (parts != cone) {
	    Domain_Free(cone);
	    value_clear(tmp);
	    Vector_Free(average);
	}
    } catch (...) {
	Domain_Free(parts);
	if (parts != cone) {
	    Domain_Free(cone);
	    value_clear(tmp);
	    Vector_Free(average);
	}
	throw;
    }
}

void barvinok_decompose(Polyhedron *C, signed_cone_consumer& scc,
			barvinok_options *options)
{
    if (options->primal)
	primal_decompose(C, scc, options);
    else
	polar_decompose(C, scc, options);
}

void vertex_decomposer::decompose_at_vertex(Param_Vertices *V, int _i, 
					    barvinok_options *options)
{
    Polyhedron *C = supporting_cone_p(P, V);
    vert = _i;
    this->V = V;

    barvinok_decompose(C, scc, options);
}

struct posneg_collector : public signed_cone_consumer {
    posneg_collector(Polyhedron *pos, Polyhedron *neg) : pos(pos), neg(neg) {}
    virtual void handle(const signed_cone& sc, barvinok_options *options) {
	Polyhedron *p = Polyhedron_Copy(sc.C);
	if (sc.sign > 0) {
	    p->next = pos;
	    pos = p;
	} else {
	    p->next = neg;
	    neg = p;
	}
    }
    Polyhedron *pos;
    Polyhedron *neg;
};

/*
 * Barvinok's Decomposition of a simplicial cone
 *
 * Returns two lists of polyhedra
 */
void barvinok_decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg)
{
    barvinok_options *options = barvinok_options_new_with_defaults();
    posneg_collector pc(*ppos, *pneg);
    POL_ENSURE_VERTICES(C);
    mat_ZZ r;
    rays(C, r);
    decompose(signed_cone(C, r, 1), pc, options);
    *ppos = pc.pos;
    *pneg = pc.neg;
    barvinok_options_free(options);
}

