#include <assert.h>
#include <math.h>
#include <glpk.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "polysign.h"

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define EMPTY_DOMAIN	-2

static LPX *solve_lp(int dir, Matrix *C, Value *f, Value denom)
{
    LPX *lp;
    int *ind;
    double *val;
    int j, k, l;
    unsigned dim = C->NbColumns-2;

    ind = ALLOCN(int, 1+dim);
    val = ALLOCN(double, 1+dim);
    lp = lpx_create_prob();
    lpx_set_obj_dir(lp, dir);
    lpx_add_rows(lp, C->NbRows);
    lpx_add_cols(lp, dim);

    for (j = 0; j < C->NbRows; ++j) {
	int type = value_zero_p(C->p[j][0]) ? LPX_FX : LPX_LO;
	for (k = 0, l = 0; k < dim; ++k) {
	    if (value_zero_p(C->p[j][1+k]))
		continue;
	    ind[1+l] = 1+k;
	    val[1+l] = VALUE_TO_DOUBLE(C->p[j][1+k]);
	    ++l;
	}
	lpx_set_mat_row(lp, 1+j, l, ind, val);
	lpx_set_row_bnds(lp, 1+j, type,
			 -VALUE_TO_DOUBLE(C->p[j][1+dim]), 0);
    }
    for (k = 0, l = 0; k < dim; ++k) {
	lpx_set_col_bnds(lp, 1+k, LPX_FR, 0, 0);
    }
    free(ind);
    free(val);
    lpx_set_int_parm(lp, LPX_K_MSGLEV, 0);

    /* objective function */
    for (j = 0; j < dim; ++j)
	lpx_set_obj_coef(lp, 1+j, VALUE_TO_DOUBLE(f[j]) /
				    VALUE_TO_DOUBLE(denom));
    lpx_set_obj_coef(lp, 0, VALUE_TO_DOUBLE(f[dim]) /
				VALUE_TO_DOUBLE(denom));

    lpx_adv_basis(lp);
    lpx_simplex(lp);

    return lp;
}

static enum lp_result constraints_affine_minmax(int dir, Matrix *C,
					      Value *f, Value denom, Value *opt)
{
    enum lp_result res = lp_ok;
    LPX *lp = solve_lp(dir, C, f, denom);

    switch (lpx_get_status(lp)) {
    case LPX_OPT:
	if (dir == LPX_MIN)
	    value_set_si(*opt, (int)ceil(lpx_get_obj_val(lp)-1e-10));
	else
	    value_set_si(*opt, (int)floor(lpx_get_obj_val(lp)+1e-10));
	break;
    case LPX_UNBND:
	res = lp_unbounded;
	break;
    case LPX_NOFEAS:
	res = lp_empty;
	break;
    default:
	assert(0);
    }
    lpx_delete_prob(lp);
    return res;
}

static int constraints_affine_minmax_sign(int dir, Matrix *C, Matrix *T,
					 int rational)
{
    LPX *lp;
    int sign;
    double opt;
    unsigned dim = C->NbColumns-2;
    assert(dim == T->NbColumns-1);
    assert(T->NbRows == 2);

    lp = solve_lp(dir, C, T->p[0], T->p[1][dim]);
    switch (lpx_get_status(lp)) {
    case LPX_OPT:
	opt = lpx_get_obj_val(lp);
	if (rational) {
	    sign = opt < 0 ? -1 : opt > 0 ? 1 : 0;
	} else {
	    if (opt < -0.5/VALUE_TO_DOUBLE(T->p[1][dim]))
		sign = -1;
	    else if (opt > 0.5/VALUE_TO_DOUBLE(T->p[1][dim]))
		sign = 1;
	    else
		sign = 0;
	}
	break;
    case LPX_UNBND:
	if (dir == LPX_MIN)
	    sign = -1;
	else
	    sign = 1;
	break;
    case LPX_NOFEAS:
	sign = EMPTY_DOMAIN;
	break;
    default:
	assert(0);
    }
    lpx_delete_prob(lp);
    return sign;
}

enum order_sign glpk_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    int rational = !POL_ISSET(options->MaxRays, POL_INTEGER);
    Matrix M;

    if (emptyQ2(D))
	return order_undefined;

    Polyhedron_Matrix_View(D, &M, D->NbConstraints);
    int min = constraints_affine_minmax_sign(LPX_MIN, &M, T, rational);
    if (min == EMPTY_DOMAIN)
	return order_undefined;
    if (min > 0)
	return order_gt;
    int max = constraints_affine_minmax_sign(LPX_MAX, &M, T, rational);
    assert(max != EMPTY_DOMAIN);
    if (max < 0)
	return order_lt;
    if (min == max)
	return order_eq;
    if (max == 0)
	return order_le;
    if (min == 0)
	return order_ge;
    return order_unknown;
}

enum lp_result glpk_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    enum lp_result res;
    Matrix M;

    if (emptyQ2(D))
	return lp_empty;

    Polyhedron_Matrix_View(D, &M, D->NbConstraints);
    res = constraints_affine_minmax(LPX_MIN, &M, obj, denom, min);
    if (res != lp_ok)
	return res;
    res = constraints_affine_minmax(LPX_MAX, &M, obj, denom, max);
    return res;
}
