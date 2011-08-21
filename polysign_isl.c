#include <assert.h>
#include <isl/lp.h>
#include <isl_set_polylib.h>
#include "polysign.h"

static void values2isl(Value *src, isl_int *dst, int len)
{
	int i;

	for (i = 0; i < len; ++i)
		isl_int_set_gmp(dst[i], src[i]);
}

static __isl_give isl_mat *extract_equalities(isl_ctx *ctx, Matrix *M)
{
	int i, j;
	int n;
	isl_int v;
	isl_mat *eq;

	n = 0;
	for (i = 0; i < M->NbRows; ++i)
		if (value_zero_p(M->p[i][0]))
			++n;

	isl_int_init(v);
	eq = isl_mat_alloc(ctx, n, M->NbColumns - 1);
	for (i = 0; i < M->NbRows; ++i) {
		if (!value_zero_p(M->p[i][0]))
			continue;
		for (j = 0; j < M->NbColumns - 1; ++j) {
			isl_int_set_gmp(v, M->p[i][1 + j]);
			eq = isl_mat_set_element(eq, i, j, v);
		}
	}
	isl_int_clear(v);

	return eq;
}

static __isl_give isl_mat *extract_inequalities(isl_ctx *ctx, Matrix *M)
{
	int i, j;
	int n;
	isl_int v;
	isl_mat *ineq;

	n = 0;
	for (i = 0; i < M->NbRows; ++i)
		if (!value_zero_p(M->p[i][0]))
			++n;

	isl_int_init(v);
	ineq = isl_mat_alloc(ctx, n, M->NbColumns - 1);
	for (i = 0; i < M->NbRows; ++i) {
		if (value_zero_p(M->p[i][0]))
			continue;
		for (j = 0; j < M->NbColumns - 1; ++j) {
			isl_int_set_gmp(v, M->p[i][1 + j]);
			ineq = isl_mat_set_element(ineq, i, j, v);
		}
	}
	isl_int_clear(v);

	return ineq;
}

enum order_sign isl_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
	isl_ctx *ctx = isl_ctx_alloc();
	isl_space *dim;
	isl_vec *aff;
	isl_basic_set *bset;
	isl_int denom;
	isl_int min, max;
	enum isl_lp_result lp_min, lp_max;
	enum order_sign sign = order_undefined;

	isl_int_init(denom);
	isl_int_init(min);
	isl_int_init(max);

	assert(D->Dimension == T->NbColumns - 1);

	aff = isl_vec_alloc(ctx, T->NbColumns);
	assert(aff);
	values2isl(T->p[0], aff->el + 1, T->NbColumns - 1);
	values2isl(T->p[0] + T->NbColumns - 1, aff->el, 1);
	values2isl(T->p[1] + T->NbColumns - 1, &denom, 1);

	dim = isl_space_set_alloc(ctx, 0, D->Dimension);
	bset = isl_basic_set_new_from_polylib(D, dim);

	lp_min = isl_basic_set_solve_lp(bset, 0, aff->el, denom, &min,
					NULL, NULL);
	assert(lp_min != isl_lp_error);

	if (lp_min == isl_lp_empty)
		sign = order_undefined;
	else if (lp_min == isl_lp_ok && isl_int_is_pos(min))
		sign = order_gt;
	else {
		lp_max = isl_basic_set_solve_lp(bset, 1, aff->el, denom, &max,
						NULL, NULL);
		assert(lp_max != isl_lp_error);

		if (lp_max == isl_lp_ok && isl_int_is_neg(max))
			sign = order_lt;
		else if (lp_min == isl_lp_ok && lp_max == isl_lp_ok &&
			 isl_int_is_zero(min) && isl_int_is_zero(max))
			sign = order_eq;
		else if (lp_min == isl_lp_ok && isl_int_is_zero(min))
			sign = order_ge;
		else if (lp_max == isl_lp_ok && isl_int_is_zero(max))
			sign = order_le;
		else
			sign = order_unknown;
	}

	isl_basic_set_free(bset);
	isl_vec_free(aff);
	isl_int_clear(min);
	isl_int_clear(max);
	isl_int_clear(denom);
	isl_ctx_free(ctx);

	return sign;
}

static enum lp_result isl_lp_result2lp_result(enum isl_lp_result res)
{
	switch (res) {
	case isl_lp_ok:
		return lp_ok;
	case isl_lp_unbounded:
		return lp_unbounded;
	case isl_lp_empty:
		return lp_empty;
	default:
		assert(0);
	}
}

enum lp_result isl_constraints_opt(Matrix *C, Value *obj, Value denom,
				    enum lp_dir dir, Value *opt)
{
	isl_ctx *ctx = isl_ctx_alloc();
	isl_space *dim;
	isl_mat *eq, *ineq;
	isl_basic_set *bset;
	isl_int isl_denom, isl_opt;
	isl_vec *aff;
	enum isl_lp_result res;
	int max = dir == lp_max;

	isl_int_init(isl_denom);
	isl_int_init(isl_opt);

	eq = extract_equalities(ctx, C);
	ineq = extract_inequalities(ctx, C);
	dim = isl_space_set_alloc(ctx, 0, C->NbColumns - 2);
	bset = isl_basic_set_from_constraint_matrices(dim, eq, ineq,
			isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
	aff = isl_vec_alloc(ctx, C->NbColumns - 1);
	assert(aff);
	values2isl(obj, aff->el + 1, aff->size - 1);
	values2isl(obj + aff->size - 1, aff->el, 1);
	isl_int_set_gmp(isl_denom, denom);

	res = isl_basic_set_solve_lp(bset, max, aff->el, isl_denom, &isl_opt,
					NULL, NULL);
	isl_int_get_gmp(isl_opt, *opt);

	isl_vec_free(aff);
	isl_int_clear(isl_denom);
	isl_int_clear(isl_opt);
	isl_basic_set_free(bset);
	isl_ctx_free(ctx);

	return isl_lp_result2lp_result(res);
}
