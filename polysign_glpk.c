#include <math.h>
#include <glpk.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "polysign.h"

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define EMPTY_DOMAIN	-2

static LPX *solve_lp(int dir, Polyhedron *P, Value *f, Value denom)
{
    LPX *lp;
    int *ind;
    double *val;
    int j, k, l;

    ind = ALLOCN(int, 1+P->Dimension);
    val = ALLOCN(double, 1+P->Dimension);
    lp = lpx_create_prob();
    lpx_set_obj_dir(lp, dir);
    lpx_add_rows(lp, P->NbConstraints);
    lpx_add_cols(lp, P->Dimension);

    for (j = 0; j < P->NbConstraints; ++j) {
	int type = j < P->NbEq ? LPX_FX : LPX_LO;
	for (k = 0, l = 0; k < P->Dimension; ++k) {
	    if (value_zero_p(P->Constraint[j][1+k]))
		continue;
	    ind[1+l] = 1+k;
	    val[1+l] = VALUE_TO_DOUBLE(P->Constraint[j][1+k]);
	    ++l;
	}
	lpx_set_mat_row(lp, 1+j, l, ind, val);
	lpx_set_row_bnds(lp, 1+j, type,
			 -VALUE_TO_DOUBLE(P->Constraint[j][1+P->Dimension]), 0);
    }
    for (k = 0, l = 0; k < P->Dimension; ++k) {
	lpx_set_col_bnds(lp, 1+k, LPX_FR, 0, 0);
    }
    free(ind);
    free(val);
    lpx_set_int_parm(lp, LPX_K_MSGLEV, 0);

    /* objective function */
    for (j = 0; j < P->Dimension; ++j)
	lpx_set_obj_coef(lp, 1+j, VALUE_TO_DOUBLE(f[j]) /
				    VALUE_TO_DOUBLE(denom));
    lpx_set_obj_coef(lp, 0, VALUE_TO_DOUBLE(f[P->Dimension]) /
				VALUE_TO_DOUBLE(denom));

    lpx_adv_basis(lp);
    lpx_simplex(lp);

    return lp;
}

static enum lp_result polyhedron_affine_minmax(int dir, Polyhedron *P,
					      Value *f, Value denom, Value *opt)
{
    enum lp_result res = lp_ok;
    LPX *lp = solve_lp(dir, P, f, denom);

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

static int polyhedron_affine_minmax_sign(int dir, Polyhedron *P, Matrix *T,
					 int rational)
{
    LPX *lp;
    int sign;
    double opt;
    assert(P->Dimension == T->NbColumns-1);
    assert(T->NbRows == 2);

    lp = solve_lp(dir, P, T->p[0], T->p[1][P->Dimension]);
    switch (lpx_get_status(lp)) {
    case LPX_OPT:
	opt = lpx_get_obj_val(lp);
	if (rational) {
	    sign = opt < 0 ? -1 : opt > 0 ? 1 : 0;
	} else {
	    if (opt < -0.5/VALUE_TO_DOUBLE(T->p[1][P->Dimension]))
		sign = -1;
	    else if (opt > 0.5/VALUE_TO_DOUBLE(T->p[1][P->Dimension]))
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

    if (emptyQ2(D))
	return order_undefined;

    int min = polyhedron_affine_minmax_sign(LPX_MIN, D, T, rational);
    if (min == EMPTY_DOMAIN)
	return order_undefined;
    if (min > 0)
	return order_gt;
    int max = polyhedron_affine_minmax_sign(LPX_MAX, D, T, rational);
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

    if (emptyQ2(D))
	return lp_empty;

    res = polyhedron_affine_minmax(LPX_MIN, D, obj, denom, min);
    if (res != lp_ok)
	return res;
    res = polyhedron_affine_minmax(LPX_MAX, D, obj, denom, max);
    return res;
}
