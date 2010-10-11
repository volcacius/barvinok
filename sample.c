#include <assert.h>
#include <isl_set_polylib.h>
#include <barvinok/util.h>
#include <barvinok/basis_reduction.h>
#include <barvinok/sample.h>
#include <barvinok/options.h>
#include "polysign.h"

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

/* If P has no rays, then we return NULL.
 * Otherwise, look for the coordinate axis with the smallest maximal non-zero
 * coefficient over all rays and a constraint that bounds the values on
 * this axis to the maximal value over the vertices plus the above maximal
 * non-zero coefficient times the number of rays minus 1.
 * Any integer point x outside this region is the sum of a point inside
 * the region and an integer multiple of the rays.
 * Write x = \sum_i a_i v_i + \sum_j b_j r_j
 * with \sum_i a_i = 1.
 * Then x = \sum_i a_i v_i + \sum_j {b_j} r_j + \sum_j [b_j] r_j
 * with y = \sum_i a_i v_i + \sum_j {b_j} r_j a point inside the region.
 */
static Polyhedron *remove_ray(Polyhedron *P, unsigned MaxRays)
{
    int r = 0;
    Vector *min, *max, *c;
    int i;
    Value s, v, tmp;
    int pos;
    Polyhedron *R;
    int rays;

    POL_ENSURE_VERTICES(P);
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

    rays = P->NbBid;
    for (r = P->NbBid; r < P->NbRays; ++r) {
	if (value_notzero_p(P->Ray[r][P->Dimension+1]))
	    continue;
	for (i = 0 ; i < P->Dimension; ++i) {
	    if (value_gt(P->Ray[r][1+i], max->p[i]))
		value_assign(max->p[i], P->Ray[r][1+i]);
	    if (value_lt(P->Ray[r][1+i], min->p[i]))
		value_assign(min->p[i], P->Ray[r][1+i]);
	}
	++rays;
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

    value_set_si(tmp, rays);
    value_addmul(v, tmp, s);
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
	remove[i] = -1;
    NbEq = 0;
    for (i = 0; i < P->NbEq; ++i) {
	int pos = First_Non_Zero(P->Constraint[i]+1, dim);
	if (First_Non_Zero(P->Constraint[i]+1+pos+1, dim-pos-1) != -1)
	    continue;
	remove[pos] = i;
	++NbEq;
    }
    assert(NbEq > 0);
    Q = Polyhedron_Alloc(P->Dimension-NbEq, P->NbConstraints-NbEq, P->NbRays);
    Q->NbEq = P->NbEq - NbEq;
    for (i = 0, k = 0; i < P->NbConstraints; ++i) {
	if (i < P->NbEq) {
	    int pos = First_Non_Zero(P->Constraint[i]+1, dim);
	    if (First_Non_Zero(P->Constraint[i]+1+pos+1, dim-pos-1) == -1)
		continue;
	}
	value_assign(Q->Constraint[k][0], P->Constraint[i][0]);
	for (j = 0, n = 0; j < P->Dimension; ++j) {
	    if (remove[j] != -1)
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
	    if (remove[j] != -1)
		++n;
	    else
		value_assign(Q->Ray[i][1+j-n], P->Ray[i][1+j]);
	}
	value_assign(Q->Ray[i][1+j-n], P->Ray[i][1+j]);
    }
    *T = Matrix_Alloc(P->Dimension+1, Q->Dimension+1);
    for (i = 0, n = 0; i < P->Dimension; ++i) {
	if (remove[i] != -1) {
	    value_oppose((*T)->p[i][Q->Dimension],
			 P->Constraint[remove[i]][1+P->Dimension]);
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

static Polyhedron *remove_all_equalities(Polyhedron *P, Matrix **T,
					 unsigned MaxRays)
{
    Matrix M;
    Polyhedron_Matrix_View(P, &M, P->NbEq);

    *T = compress_variables(&M, 0);

    if (!*T)
	return NULL;
    P = Polyhedron_Preimage(P, *T, MaxRays);

    return P;
}

static Vector *product_sample(Polyhedron *P, Matrix *T,
			      struct barvinok_options *options)
{
    int i;
    Vector *sample = NULL;
    Vector *tmp = Vector_Alloc(T->NbRows);
    i = 0;
    for (; P; P = P->next) {
	Vector *P_sample;
	Polyhedron *next = P->next;
	P->next = NULL;
	P_sample = Polyhedron_Sample(P, options);
	P->next = next;
	if (!P_sample) {
	    Vector_Free(tmp);
	    tmp = NULL;
	    break;
	}
	Vector_Copy(P_sample->p, tmp->p+i, P->Dimension);
	Vector_Free(P_sample);
	i += P->Dimension;
    }
    if (tmp) {
	sample = Vector_Alloc(T->NbRows + 1);
	Matrix_Vector_Product(T, tmp->p, sample->p);
	value_set_si(sample->p[T->NbRows], 1);
	Vector_Free(tmp);
    }
    return sample;
}

static Vector *isl_Polyhedron_Sample(Polyhedron *P,
	struct barvinok_options *options)
{
	int i;
	isl_ctx *ctx = isl_ctx_alloc();
	isl_dim *dim;
	int nvar = P->Dimension;
	isl_basic_set *bset;
	isl_point *pnt;
	Vector *sample = NULL;

	dim = isl_dim_set_alloc(ctx, 0, nvar);
	bset = isl_basic_set_new_from_polylib(P, dim);
	pnt = isl_basic_set_sample_point(bset);

	if (!isl_point_is_void(pnt)) {
		isl_int v;

		isl_int_init(v);
		sample = Vector_Alloc(1 + nvar);
		assert(sample);
		for (i = 0; i < nvar; ++i) {
			isl_point_get_coordinate(pnt, isl_dim_set, i, &v);
			isl_int_get_gmp(v, sample->p[i]);
		}
		value_set_si(sample->p[nvar], 1);
		isl_int_clear(v);
	}

	isl_point_free(pnt);

	isl_ctx_free(ctx);

	return sample;
}

/* This function implements the algorithm described in
 * "An Implementation of the Generalized Basis Reduction Algorithm
 *  for Integer Programming" of Cook el al. to find an integer point
 * in a polyhedron.
 * If the polyhedron is unbounded, we first remove its rays.
 */
Vector *Polyhedron_Sample(Polyhedron *P, struct barvinok_options *options)
{
    int i, j;
    Vector *sample = NULL, *obj = NULL;
    Polyhedron *Q;
    Matrix *inv = NULL;
    Value min, max, tmp;
    Vector *v;
    int ok;
    enum lp_result res;
    Matrix *T;

    if (options->gbr_lp_solver == BV_GBR_ISL)
	return isl_Polyhedron_Sample(P, options);

    if (emptyQ2(P))
	return NULL;

    if (P->Dimension == 0) {
	sample = Vector_Alloc(1);
	value_set_si(sample->p[0], 1);
	return sample;
    }

    if (P->Dimension == 1)
	POL_ENSURE_VERTICES(P);

    for (i = 0; i < P->NbRays; ++i)
	if (value_one_p(P->Ray[i][1+P->Dimension])) {
	    sample = Vector_Alloc(P->Dimension+1);
	    Vector_Copy(P->Ray[i]+1, sample->p, P->Dimension+1);
	    return sample;
	}

    if (P->NbEq > 0) {
	Vector *Q_sample;
	Polyhedron *Q = remove_all_equalities(P, &T, options->MaxRays);
	if (!Q)
	    return NULL;
	Q_sample = Polyhedron_Sample(Q, options);
	Polyhedron_Free(Q);
	if (Q_sample) {
	    sample = Vector_Alloc(P->Dimension + 1);
	    Matrix_Vector_Product(T, Q_sample->p, sample->p);
	    Vector_Free(Q_sample);
	}
	Matrix_Free(T);
	return sample;
    }

    Q = Polyhedron_Factor(P, 0, &T, options->MaxRays);
    if (Q) {
	sample = product_sample(Q, T, options);
	Domain_Free(Q);
	Matrix_Free(T);
	return sample;
    }

    value_init(min);
    value_init(max);

    obj = Vector_Alloc(P->Dimension+1);
    value_set_si(obj->p[0], 1);
    res = polyhedron_range(P, obj->p, obj->p[0], &min, &max, options);
    if (res == lp_unbounded) {
unbounded:
	value_clear(min);
	value_clear(max);
	Vector_Free(obj);

	Q = remove_ray(P, options->MaxRays);
	assert(Q);
	sample = Polyhedron_Sample(Q, options);
	Polyhedron_Free(Q);
	return sample;
    }
    if (res == lp_empty)
	goto out;
    assert(res == lp_ok);

    if (value_eq(min, max)) {
	Q = P;
    } else {
	Matrix *M, *T;
	Matrix *basis;

	options->gbr_only_first = 1;
	basis = Polyhedron_Reduced_Basis(P, options);
	options->gbr_only_first = 0;

	if (!basis)
	    goto unbounded;
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

	Q = Polyhedron_Image(P, T, options->MaxRays);
	res = polyhedron_range(Q, obj->p, obj->p[0], &min, &max, options);

	Matrix_Free(T);
	if (res == lp_empty)
	    goto out;
	assert(res == lp_ok);
    }

    value_init(tmp);

    v = Vector_Alloc(1+Q->Dimension+1);
    value_set_si(v->p[1], -1);

    for (value_assign(tmp, min); value_le(tmp, max); value_increment(tmp, tmp)) {
	Polyhedron *R, *S;
	Matrix *T;
	Vector *S_sample;
	value_assign(v->p[1+Q->Dimension], tmp);

	R = AddConstraints(v->p, 1, Q, options->MaxRays);
	R = DomainConstraintSimplify(R, options->MaxRays);
	if (emptyQ(R)) {
	    Polyhedron_Free(R);
	    continue;
	}

	S = Polyhedron_RemoveFixedColumns(R, &T);
	Polyhedron_Free(R);
	S_sample = Polyhedron_Sample(S, options);
	Polyhedron_Free(S);
	if (S_sample) {
	    Vector *Q_sample = obj;
	    obj = NULL;
	    Matrix_Vector_Product(T, S_sample->p, Q_sample->p);
	    Matrix_Free(T);
	    Vector_Free(S_sample);
	    if (!inv)
		sample = Q_sample;
	    else {
		sample = Vector_Alloc(P->Dimension + 1);
		Matrix_Vector_Product(inv, Q_sample->p, sample->p);
		Vector_Free(Q_sample);
	    }
	    break;
	}
	Matrix_Free(T);
    }
    value_clear(tmp);

    Vector_Free(v);

out:
    if (inv)
	Matrix_Free(inv);
    if (P != Q)
	Polyhedron_Free(Q);
    if (obj)
	Vector_Free(obj);
    value_clear(min);
    value_clear(max);

    return sample;
}
