#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include "basis_reduction.h"
#include "sample.h"

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

/* If P has no rays, then we return NULL.
 * Otherwise, look for the coordinate axis with the smallest maximal non-zero
 * coefficient over all rays and a constraint that bounds the values on
 * this axis to the maximal value over the vertices plus the above maximal
 * non-zero coefficient minus 1.
 * Any integer point outside this region should be the sum of a point inside
 * the region and an integer multiple of the rays.
 */
static Polyhedron *remove_ray(Polyhedron *P, unsigned MaxRays)
{
    int r = 0;
    Vector *min, *max, *c;
    int i;
    Value s, v, tmp;
    int pos;
    Polyhedron *R;

    if (P->NbBid == 0)
	for (; r < P->NbRays; ++r)
	    if (value_zero_p(P->Ray[r][P->Dimension+1]))
		break;
    if (P->NbBid == 0 && r == P->NbRays)
	return NULL;

    max = Vector_Alloc(P->Dimension);
    min = Vector_Alloc(P->Dimension);
    for (r = 0; r < P->NbBid; ++r)
	for (i = 0 ; i < P->Dimension; ++i)
	    if (value_abs_gt(P->Ray[r][1+i], max->p[i]))
		value_absolute(max->p[i], P->Ray[r][1+i]);

    for (i = 0 ; i < P->Dimension; ++i)
	value_oppose(min->p[i], max->p[i]);

    for (r = P->NbBid; r < P->NbRays; ++r) {
	if (value_notzero_p(P->Ray[r][P->Dimension+1]))
	    continue;
	for (i = 0 ; i < P->Dimension; ++i) {
	    if (value_gt(P->Ray[r][1+i], max->p[i]))
		value_assign(max->p[i], P->Ray[r][1+i]);
	    if (value_lt(P->Ray[r][1+i], min->p[i]))
		value_assign(min->p[i], P->Ray[r][1+i]);
	}
    }

    value_init(s);
    value_init(v);
    value_init(tmp);

    for (i = 0 ; i < P->Dimension; ++i) {
	if (value_notzero_p(min->p[i]) && 
	    (value_zero_p(s) || value_abs_lt(min->p[i], s))) {
	    value_assign(s, min->p[i]);
	    pos = i;
	}
	if (value_notzero_p(max->p[i]) && 
	    (value_zero_p(s) || value_abs_lt(max->p[i], s))) {
	    value_assign(s, max->p[i]);
	    pos = i;
	}
    }

    for (r = P->NbBid; r < P->NbRays; ++r)
	if (value_notzero_p(P->Ray[r][P->Dimension+1]))
	    break;

    if (value_pos_p(s))
	mpz_cdiv_q(v, P->Ray[r][1+pos], P->Ray[r][P->Dimension+1]);
    else
	mpz_fdiv_q(v, P->Ray[r][1+pos], P->Ray[r][P->Dimension+1]);

    for ( ; r < P->NbRays; ++r) {
	if (value_zero_p(P->Ray[r][P->Dimension+1]))
	    continue;

	if (value_pos_p(s)) {
	    mpz_cdiv_q(tmp, P->Ray[r][1+pos], P->Ray[r][P->Dimension+1]);
	    if (value_gt(tmp, v))
		value_assign(v, tmp);
	} else {
	    mpz_fdiv_q(tmp, P->Ray[r][1+pos], P->Ray[r][P->Dimension+1]);
	    if (value_lt(tmp, v))
		value_assign(v, tmp);
	}
    }

    c = Vector_Alloc(1+P->Dimension+1);

    value_addto(v, v, s);
    value_set_si(c->p[0], 1);
    if (value_pos_p(s)) {
	value_set_si(c->p[1+pos], -1);
	value_assign(c->p[1+P->Dimension], v);
    } else {
	value_set_si(c->p[1+pos], 1);
	value_oppose(c->p[1+P->Dimension], v);
    }
    value_decrement(c->p[1+P->Dimension], c->p[1+P->Dimension]);

    R = AddConstraints(c->p, 1, P, MaxRays);

    Vector_Free(c);

    Vector_Free(min);
    Vector_Free(max);

    value_clear(tmp);
    value_clear(s);
    value_clear(v);

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
    int i, j, k, n;
    int dim = P->Dimension;
    int *remove = ALLOCN(int, dim);
    Polyhedron *Q;
    int NbEq;

    assert(POL_HAS(P, POL_INEQUALITIES));
    for (i = 0; i < dim; ++i)
	remove[i] = 0;
    NbEq = 0;
    for (i = 0; i < P->NbEq; ++i) {
	int pos = First_Non_Zero(P->Constraint[i]+1, dim);
	if (First_Non_Zero(P->Constraint[i]+1+pos+1, dim-pos-1) != -1)
	    continue;
	remove[pos] = 1;
	++NbEq;
    }
    assert(NbEq > 0);
    Q = Polyhedron_Alloc(P->Dimension-NbEq, P->NbConstraints-NbEq, P->NbRays);
    for (i = 0, k = 0; i < P->NbConstraints; ++i) {
	if (i < P->NbEq) {
	    int pos = First_Non_Zero(P->Constraint[i]+1, dim);
	    if (First_Non_Zero(P->Constraint[i]+1+pos+1, dim-pos-1) == -1)
		continue;
	}
	value_assign(Q->Constraint[k][0], P->Constraint[i][0]);
	for (j = 0, n = 0; j < P->Dimension; ++j) {
	    if (remove[j])
		++n;
	    else
		value_assign(Q->Constraint[k][1+j-n], P->Constraint[i][1+j]);
	}
	value_assign(Q->Constraint[k][1+j-n], P->Constraint[i][1+j]);
	++k;
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
    if (emptyQ(P))
	return NULL;

    if (P->Dimension == 0) {
	sample = Vector_Alloc(1);
	value_set_si(sample->p[0], 1);
	return sample;
    }

    for (i = 0; i < P->NbRays; ++i)
	if (value_one_p(P->Ray[i][1+P->Dimension])) {
	    sample = Vector_Alloc(P->Dimension+1);
	    Vector_Copy(P->Ray[i]+1, sample->p, P->Dimension+1);
	    return sample;
	}

    Q = remove_ray(P, MaxRays);
    if (Q) {
	sample = Polyhedron_Sample(Q, MaxRays);
	Polyhedron_Free(Q);
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

    POL_ENSURE_VERTICES(Q);

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
	R = DomainConstraintSimplify(R, MaxRays);
	if (emptyQ(R)) {
	    Polyhedron_Free(R);
	    continue;
	}

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
