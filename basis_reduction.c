#include <math.h>
#include <glpk.h>
#include <polylib/polylibgmp.h>

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static LPX *init_lp(Polyhedron *P)
{
    LPX *lp;
    int *ind;
    double *val;
    int i, j, k, l;

    ind = ALLOCN(int, 1+P->Dimension);
    val = ALLOCN(double, 1+P->Dimension);
    lp = lpx_create_prob();
    lpx_set_obj_dir(lp, LPX_MAX);
    lpx_add_rows(lp, 2*P->NbConstraints);
    lpx_add_cols(lp, 2*P->Dimension);
    for (i = 0; i < 2; ++i) {
	for (j = 0; j < P->NbConstraints; ++j) {
	    for (k = 0, l = 0; k < P->Dimension; ++k) {
		if (value_zero_p(P->Constraint[j][1+k]))
		    continue;
		ind[1+l] = 1+i*P->Dimension+k;
		val[1+l] = VALUE_TO_DOUBLE(P->Constraint[j][1+k]);
		++l;
	    }
	    lpx_set_mat_row(lp, 1+i*P->NbConstraints+j, l, ind, val);
	    lpx_set_row_bnds(lp, 1+i*P->NbConstraints+j, LPX_LO,
			     -VALUE_TO_DOUBLE(P->Constraint[j][1+P->Dimension]), 0);
	}
	for (k = 0, l = 0; k < P->Dimension; ++k) {
	    lpx_set_col_bnds(lp, 1+i*P->Dimension+k, LPX_FR, 0, 0);
	}
    }
    free(ind);
    free(val);
    return lp;
}

static void set_lp_obj(LPX *lp, Value *row, int dim)
{
    int j;
    for (j = 0; j < dim; ++j) {
	lpx_set_obj_coef(lp, 1+j, VALUE_TO_DOUBLE(row[j]));
	lpx_set_obj_coef(lp, 1+dim+j, -VALUE_TO_DOUBLE(row[j]));
    }
}

static int add_lp_row(LPX *lp, Value *row, int dim)
{
    int j, l;
    int nr = lpx_add_rows(lp, 1);
    int *ind;
    double *val;

    ind = ALLOCN(int, 1+2*dim);
    val = ALLOCN(double, 1+2*dim);
    for (j = 0, l = 0; j < dim; ++j) {
	if (value_zero_p(row[j]))
	    continue;
	ind[1+l] = 1+j;
	val[1+l] = VALUE_TO_DOUBLE(row[j]);
	ind[1+l+1] = 1+dim+j;
	val[1+l+1] = -VALUE_TO_DOUBLE(row[j]);
	l += 2;
    }
    lpx_set_mat_row(lp, nr, l, ind, val);
    lpx_set_row_bnds(lp, nr, LPX_FX, 0, 0);
    free(ind);
    free(val);

    return nr;
}

static void del_lp_row(LPX *lp)
{
    int rows[2];
    rows[1] = lpx_get_num_rows(lp);
    lpx_del_rows(lp, 1, rows);
}

static void save_alpha(LPX *lp, int first, int n, double *alpha)
{
    int i;

    for (i = 0; i < n; ++i)
	alpha[i] = -lpx_get_row_dual(lp, first+i);
}

/* This function implements the algorithm described in
 * "An Implementation of the Generalized Basis Reduction Algorithm
 *  for Integer Programming" of Cook el al. to compute a reduced basis.
 * We use \epsilon = 1/4.
 */
Matrix *reduced_basis(Polyhedron *P)
{
    int dim = P->Dimension;
    int i;
    Matrix *basis = Identity(dim);
    LPX *lp;
    double F_old, alpha, F_new;
    int row;
    Value one, tmp;
    Vector *b_tmp;
    double *F;
    double *alpha_buffer[2];
    double *alpha_saved;
    double F_saved;
    int use_saved = 0;

    if (P->Dimension == 1)
	return basis;

    value_init(one);
    value_init(tmp);
    value_set_si(one, 1);

    b_tmp = Vector_Alloc(dim);

    F = ALLOCN(double, dim);
    alpha_buffer[0] = ALLOCN(double, dim);
    alpha_buffer[1] = ALLOCN(double, dim);
    alpha_saved = alpha_buffer[0];

    lp = init_lp(P);

    lpx_set_int_parm(lp, LPX_K_MSGLEV, 0);

    i = 0;

    set_lp_obj(lp, basis->p[0], dim);

    lpx_adv_basis(lp);
    lpx_simplex(lp);
    F[0] = lpx_get_obj_val(lp);
    assert(F[0] > -1e-10);
    if (F[0] < 0)
	F[0] = 0;

    do {
	int mu[2];

	if (use_saved) {
	    row = lpx_get_num_rows(lp)+1;
	    F_new = F_saved;
	    alpha = alpha_saved[i];
	} else {
	    row = add_lp_row(lp, basis->p[i], dim);
	    set_lp_obj(lp, basis->p[i+1], dim);
	    lpx_adv_basis(lp);
	    lpx_simplex(lp);
	    F_new = lpx_get_obj_val(lp);

	    alpha = -lpx_get_row_dual(lp, row);

	    if (i > 0)
		save_alpha(lp, row-i, i, alpha_saved);

	    del_lp_row(lp);
	}
	F[i+1] = F_new;
	assert(F[i+1] > -1e-10);
	if (F[i+1] < 0)
	    F[i+1] = 0;


	mu[0] = (int)floor(alpha+1e-10);
	mu[1] = (int)ceil(alpha-1e-10);

	if (mu[0] == mu[1])
	    value_set_si(tmp, mu[0]);
	else {
	    double mu_F[2];
	    int j;

	    for (j = 0; j <= 1; ++j) {
		value_set_si(tmp, mu[j]);
		Vector_Combine(basis->p[i+1], basis->p[i], b_tmp->p, one, tmp, dim);
		set_lp_obj(lp, b_tmp->p, dim);
		lpx_adv_basis(lp);
		lpx_simplex(lp);
		mu_F[j] = lpx_get_obj_val(lp);
		if (i > 0)
		    save_alpha(lp, row-i, i, alpha_buffer[j]);
	    }

	    if (mu_F[0] < mu_F[1])
		j = 0;
	    else
		j = 1;

	    value_set_si(tmp, mu[j]);
	    F_new = mu_F[j];
	    alpha_saved = alpha_buffer[j];
	}
	Vector_Combine(basis->p[i+1], basis->p[i], basis->p[i+1], one, tmp, dim);

	assert(F_new > -1e-10);
	if (F_new < 0)
	    F_new = 0;
	F_old = F[i];

	use_saved = 0;
	if (F_new < 3*F_old/4) {
	    Vector_Exchange(basis->p[i], basis->p[i+1], dim);
	    if (i > 0) {
		use_saved = 1;
		F_saved = F_new;
		del_lp_row(lp);
		--i;
	    } else
		F[0] = F_new;
	} else {
	    add_lp_row(lp, basis->p[i], dim);
	    ++i;
	}
    } while (i < dim-1);

    Vector_Free(b_tmp);

    value_clear(one);
    value_clear(tmp);
    free(F);
    free(alpha_buffer[0]);
    free(alpha_buffer[1]);

    lpx_delete_prob(lp);

    return basis;
}
