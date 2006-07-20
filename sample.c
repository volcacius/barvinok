#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include "basis_reduction.h"
#include "sample.h"

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

/* If P has no rays, then we return NULL.
 * Otherwise, transform the polyhedron such that one of the rays
 * is the first unit vector and cut it off at a height that ensures
 * that if the whole polyhedron has any points, then the remaining part
 * has integer points.  In particular we add the largest coefficient
 * of a ray to the highest vertex (rounded up).
 *
 * The matrix that transforms the resulting polytope to part of the
 * original polyhedron is returned through T.
 */
static Polyhedron *remove_ray(Polyhedron *P, Matrix **T, unsigned MaxRays)
{
    int r = 0;
    Matrix *M, *M2;
    Value c, tmp;
    Value g;
    int first;
    Vector *v;
    Value offset, size;
    Polyhedron *R;
    int i;

    if (P->NbBid == 0)
	for (; r < P->NbRays; ++r)
	    if (value_zero_p(P->Ray[r][P->Dimension+1]))
		break;
    if (P->NbBid == 0 && r == P->NbRays)
	return NULL;

    value_init(g);
    v = Vector_Alloc(P->Dimension+1);
    Vector_Gcd(P->Ray[r]+1, P->Dimension, &g);
    Vector_AntiScale(P->Ray[r]+1, v->p, g, P->Dimension+1);
    M = unimodular_complete(v);
    value_set_si(M->p[P->Dimension][P->Dimension], 1);
    M2 = Transpose(M);
    Matrix_Free(M);
    P = Polyhedron_Preimage(P, M2, 0);
    *T = M2;
    value_clear(g);
    Vector_Free(v);

    first = 1;
    value_init(offset);
    value_init(size);
    value_init(tmp);
    value_set_si(size, 0);

    for (i = 0; i < P->NbBid; ++i) {
	value_absolute(tmp, P->Ray[i][1]);
	if (value_gt(tmp, size))
	    value_assign(size, tmp);
    }
    for (i = P->NbBid; i < P->NbRays; ++i) {
	if (value_zero_p(P->Ray[i][P->Dimension+1])) {
	    if (value_gt(P->Ray[i][1], size))
		value_assign(size, P->Ray[i][1]);
	    continue;
	}
	mpz_cdiv_q(tmp, P->Ray[i][1], P->Ray[i][P->Dimension+1]);
	if (first || value_gt(tmp, offset)) {
	    value_assign(offset, tmp);
	    first = 0;
	}
    }
    value_addto(offset, offset, size);
    value_clear(size);
    value_clear(tmp);

    v = Vector_Alloc(P->Dimension+2);
    value_set_si(v->p[0], 1);
    value_set_si(v->p[1], -1);
    value_assign(v->p[1+P->Dimension], offset);
    R = AddConstraints(v->p, 1, P, MaxRays);
    value_clear(offset);
    Vector_Free(v);
    Polyhedron_Free(P);

    return R;
}

static void print_minmax(Polyhedron *P)
{
    int i, j;
    POL_ENSURE_VERTICES(P);
    Polyhedron_Print(stderr, P_VALUE_FMT, P);
    for (i = 0; i < P->Dimension; ++i) {
	Value min, max, tmp;
	value_init(min);
	value_init(max);
	value_init(tmp);

	mpz_cdiv_q(min, P->Ray[0][1+i], P->Ray[0][1+P->Dimension]);
	mpz_fdiv_q(max, P->Ray[0][1+i], P->Ray[0][1+P->Dimension]);

	for (j = 1; j < P->NbRays; ++j) {
	    mpz_cdiv_q(tmp, P->Ray[j][1+i], P->Ray[j][1+P->Dimension]);
	    if (value_lt(tmp, min))
		value_assign(min, tmp);
	    mpz_fdiv_q(tmp, P->Ray[j][1+i], P->Ray[j][1+P->Dimension]);
	    if (value_gt(tmp, max))
		value_assign(max, tmp);
	}
	fprintf(stderr, "i: %d,  min: ", i);
	value_print(stderr, VALUE_FMT, min);
	fprintf(stderr, ", max: ");
	value_print(stderr, VALUE_FMT, max);
	fprintf(stderr, "\n");

	value_clear(min);
	value_clear(max);
	value_clear(tmp);
    }
}

/* Remove coordinates that have a fixed value and return the matrix
 * that adds these fixed coordinates again through T.
 */
static Polyhedron *Polyhedron_RemoveFixedColumns(Polyhedron *P, Matrix **T)
{
    int i, j, n;
    int dim = P->Dimension;
    int *remove = ALLOCN(int, dim);
    Polyhedron *Q;

    assert(POL_HAS(P, POL_INEQUALITIES));
    for (i = 0; i < dim; ++i)
	remove[i] = 0;
    for (i = 0; i < P->NbEq; ++i) {
	int pos = First_Non_Zero(P->Constraint[i]+1, dim);
	assert(First_Non_Zero(P->Constraint[i]+1+pos+1, dim-pos-1) == -1);
	remove[pos] = 1;
    }
    Q = Polyhedron_Alloc(P->Dimension-P->NbEq, P->NbConstraints-P->NbEq, P->NbRays);
    for (i = 0; i < Q->NbConstraints; ++i) {
	value_assign(Q->Constraint[i][0], P->Constraint[P->NbEq+i][0]);
	for (j = 0, n = 0; j < P->Dimension; ++j) {
	    if (remove[j])
		++n;
	    else
		value_assign(Q->Constraint[i][1+j-n], P->Constraint[P->NbEq+i][1+j]);
	}
	value_assign(Q->Constraint[i][1+j-n], P->Constraint[P->NbEq+i][1+j]);
    }
    for (i = 0; i < Q->NbRays; ++i) {
	value_assign(Q->Ray[i][0], P->Ray[i][0]);
	for (j = 0, n = 0; j < P->Dimension; ++j) {
	    if (remove[j])
		++n;
	    else
		value_assign(Q->Ray[i][1+j-n], P->Ray[i][1+j]);
	}
	value_assign(Q->Ray[i][1+j-n], P->Ray[i][1+j]);
    }
    *T = Matrix_Alloc(P->Dimension+1, Q->Dimension+1);
    for (i = 0, n = 0; i < P->Dimension; ++i) {
	if (remove[i]) {
	    value_oppose((*T)->p[i][Q->Dimension], P->Constraint[n][1+P->Dimension]);
	    ++n;
	} else
	    value_set_si((*T)->p[i][i-n], 1);
    }
    value_set_si((*T)->p[i][i-n], 1);
    POL_SET(Q, POL_VALID);
    if (POL_HAS(P, POL_INEQUALITIES))
	POL_SET(Q, POL_INEQUALITIES);
    if (POL_HAS(P, POL_FACETS))
	POL_SET(Q, POL_FACETS);
    if (POL_HAS(P, POL_POINTS))
	POL_SET(Q, POL_POINTS);
    if (POL_HAS(P, POL_VERTICES))
	POL_SET(Q, POL_VERTICES);
    free(remove);
    return Q;
}

/* This function implements the algorithm described in
 * "An Implementation of the Generalized Basis Reduction Algorithm
 *  for Integer Programming" of Cook el al. to find an integer point
 * in a polyhedron.
 * If the polyhedron is unbounded, we first remove its rays.
 */
Vector *Polyhedron_Sample(Polyhedron *P, unsigned MaxRays)
{
    int i, j;
    Vector *sample = NULL;
    Polyhedron *Q;
    Matrix *T, *inv, *M;
    Value min, max, tmp;
    Vector *v;
    int ok;

    POL_ENSURE_VERTICES(P);

    for (i = 0; i < P->NbRays; ++i)
	if (value_one_p(P->Ray[i][1+P->Dimension])) {
	    sample = Vector_Alloc(P->Dimension+1);
	    Vector_Copy(P->Ray[i]+1, sample->p, P->Dimension+1);
	    return sample;
	}

    /* for now */
    assert(P->NbEq == 0);

    Q = remove_ray(P, &T, MaxRays);
    if (Q) {
	Vector *Q_sample;

	Q_sample = Polyhedron_Sample(Q, MaxRays);
	Polyhedron_Free(Q);

	if (Q_sample) {
	    sample = Vector_Alloc(P->Dimension + 1);
	    Matrix_Vector_Product(T, Q_sample->p, sample->p);
	    Vector_Free(Q_sample);
	}

	Matrix_Free(T);
	return sample;
    }

    Matrix *basis = reduced_basis(P);

    T = Matrix_Alloc(P->Dimension+1, P->Dimension+1);
    inv = Matrix_Alloc(P->Dimension+1, P->Dimension+1);
    for (i = 0; i < P->Dimension; ++i)
	for (j = 0; j < P->Dimension; ++j)
	    value_assign(T->p[i][j], basis->p[i][j]);
    value_set_si(T->p[P->Dimension][P->Dimension], 1);
    Matrix_Free(basis);

    M = Matrix_Copy(T);
    ok = Matrix_Inverse(M, inv);
    assert(ok);
    Matrix_Free(M);

    Q = Polyhedron_Image(P, T, MaxRays);

    value_init(min);
    value_init(max);
    value_init(tmp);

    mpz_cdiv_q(min, Q->Ray[0][1], Q->Ray[0][1+Q->Dimension]);
    mpz_fdiv_q(max, Q->Ray[0][1], Q->Ray[0][1+Q->Dimension]);

    for (j = 1; j < Q->NbRays; ++j) {
	mpz_cdiv_q(tmp, Q->Ray[j][1], Q->Ray[j][1+Q->Dimension]);
	if (value_lt(tmp, min))
	    value_assign(min, tmp);
	mpz_fdiv_q(tmp, Q->Ray[j][1], Q->Ray[j][1+Q->Dimension]);
	if (value_gt(tmp, max))
	    value_assign(max, tmp);
    }

    v = Vector_Alloc(1+Q->Dimension+1);
    value_set_si(v->p[1], -1);

    for (value_assign(tmp, min); value_le(tmp, max); value_increment(tmp, tmp)) {
	Polyhedron *R, *S;
	Matrix *T;
	Vector *S_sample;
	value_assign(v->p[1+Q->Dimension], tmp);

	R = AddConstraints(v->p, 1, Q, MaxRays);
	S = Polyhedron_RemoveFixedColumns(R, &T);
	Polyhedron_Free(R);
	S_sample = Polyhedron_Sample(S, MaxRays);
	Polyhedron_Free(S);
	if (S_sample) {
	    Vector *Q_sample = Vector_Alloc(Q->Dimension + 1);
	    Matrix_Vector_Product(T, S_sample->p, Q_sample->p);
	    Matrix_Free(T);
	    Vector_Free(S_sample);
	    sample = Vector_Alloc(P->Dimension + 1);
	    Matrix_Vector_Product(inv, Q_sample->p, sample->p);
	    Vector_Free(Q_sample);
	    break;
	}
	Matrix_Free(T);
    }

    Matrix_Free(T);
    Matrix_Free(inv);
    Polyhedron_Free(Q);
    Vector_Free(v);

    value_clear(min);
    value_clear(max);
    value_clear(tmp);

    return sample;
}
