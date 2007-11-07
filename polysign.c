#include <assert.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "polysign.h"
#include "config.h"

static int is_unit_row(Value *row, int pos, int len)
{
    if (!value_one_p(row[pos]) && !value_mone_p(row[pos]))
	return 0;
    return First_Non_Zero(row+pos+1, len-(pos+1)) == -1;
}

/* Transform the constraints of P to "standard form".
 * In particular, if P is described by
 *		A x + b(p) >= 0
 * then this function returns a matrix H = A U, A = H Q, such
 * that D x' = D Q x >= -b(p), with D a diagonal matrix with
 * positive entries.  The calling function can then construct
 * the standard form H' x' - I s + b(p) = 0, with H' the rows of H
 * that are not positive multiples of unit vectors
 * (since those correspond to D x' >= -b(p)).
 * The number of rows in H' is returned in *rows_p.
 * Optionally returns the matrix that maps the new variables
 * back to the old variables x = U x'.
 * Note that the rows of H (and P) may be reordered by this function.
 */
Matrix *standard_constraints(Polyhedron *P, unsigned nparam, int *rows_p,
			     Matrix **T)
{
    unsigned nvar = P->Dimension - nparam;
    int i, j, d;
    int rows;
    Matrix *M;
    Matrix *H, *U, *Q;

    assert(P->NbEq == 0);

    /* move constraints only involving parameters down
     * and move unit vectors (if there are any) to the right place.
     */
    for (d = 0, j = P->NbConstraints; d < j; ++d) {
	int index;
	index = First_Non_Zero(P->Constraint[d]+1, nvar);
	if (index != -1) {
	    if (index != d &&
		is_unit_row(P->Constraint[d]+1, index, nvar)) {
		Vector_Exchange(P->Constraint[d], P->Constraint[index],
				P->Dimension+2);
		if (index > d)
		    --d;
	    }
	    continue;
	}
	while (d < --j && First_Non_Zero(P->Constraint[j]+1, nvar) == -1)
	    ;
	if (d >= j)
	    break;
	Vector_Exchange(P->Constraint[d], P->Constraint[j], P->Dimension+2);
    }
    M = Matrix_Alloc(d, nvar);
    for (j = 0; j < d; ++j)
	Vector_Copy(P->Constraint[j]+1, M->p[j], nvar);

    neg_left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    Matrix_Free(Q);
    if (T)
	*T = U;
    else
	Matrix_Free(U);

    /* Rearrange rows such that top of H is lower diagonal and
     * compute the number of non (multiple of) unit-vector rows.
     */
    rows = H->NbRows-nvar;
    for (i = 0; i < H->NbColumns; ++i) {
	for (j = i; j < H->NbRows; ++j)
	    if (value_notzero_p(H->p[j][i]))
		break;
	if (j != i) {
	    Vector_Exchange(P->Constraint[i], P->Constraint[j], P->Dimension+2);
	    Vector_Exchange(H->p[i], H->p[j], H->NbColumns);
	}
	if (First_Non_Zero(H->p[i], i) != -1)
	    rows++;
    }
    *rows_p = rows;

    return H;
}

#ifndef HAVE_LIBGLPK
enum order_sign glpk_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    assert(0);
}

enum lp_result glpk_constraints_opt(Matrix *C, Value *obj, Value denom,
				    enum lp_dir dir, Value *opt)
{
    assert(0);
}

enum lp_result glpk_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    assert(0);
}
#endif

#ifndef HAVE_LIBCDDGMP
enum order_sign cdd_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    assert(0);
}

enum order_sign cddf_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    assert(0);
}

enum lp_result cdd_constraints_opt(Matrix *C, Value *obj, Value denom,
				    enum lp_dir dir, Value *opt)
{
    assert(0);
}

enum lp_result cddf_constraints_opt(Matrix *C, Value *obj, Value denom,
				    enum lp_dir dir, Value *opt)
{
    assert(0);
}

enum lp_result cdd_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    assert(0);
}

enum lp_result cddf_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    assert(0);
}
#endif

enum order_sign polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    if (options->lp_solver == BV_LP_POLYLIB)
	return PL_polyhedron_affine_sign(D, T, options);
    else if (options->lp_solver == BV_LP_GLPK)
	return glpk_polyhedron_affine_sign(D, T, options);
    else if (options->lp_solver == BV_LP_CDD)
	return cdd_polyhedron_affine_sign(D, T, options);
    else if (options->lp_solver == BV_LP_CDDF)
	return cddf_polyhedron_affine_sign(D, T, options);
    else
	assert(0);
}

/*
 * Optimize (minimize or maximize depending on dir) the affine
 * objective function obj (of length dimension+1), with denominator
 * denom over the polyhedron specified by the constraints C.
 * The result is returned in opt.
 */
enum lp_result constraints_opt(Matrix *C, Value *obj, Value denom,
				enum lp_dir dir, Value *opt,
				struct barvinok_options *options)
{
    if (options->lp_solver == BV_LP_POLYLIB)
	return PL_constraints_opt(C, obj, denom, dir, opt, options->MaxRays);
    else if (options->lp_solver == BV_LP_GLPK)
	return glpk_constraints_opt(C, obj, denom, dir, opt);
    else if (options->lp_solver == BV_LP_CDD)
	return cdd_constraints_opt(C, obj, denom, dir, opt);
    else if (options->lp_solver == BV_LP_CDDF)
	return cddf_constraints_opt(C, obj, denom, dir, opt);
    else if (options->lp_solver == BV_LP_PIP)
	return pip_constraints_opt(C, obj, denom, dir, opt);
    else
	assert(0);
}

/*
 * Optimize (minimize or maximize depending on dir) the affine
 * objective function obj (of length dimension+1), with denominator
 * denom over the polyhedron specified by P.
 * The result is returned in opt.
 */
enum lp_result polyhedron_opt(Polyhedron *P, Value *obj, Value denom,
				enum lp_dir dir, Value *opt,
				struct barvinok_options *options)
{
    if (options->lp_solver == BV_LP_POLYLIB)
	return PL_polyhedron_opt(P, obj, denom, dir, opt);
    else {
	Matrix M;
	Polyhedron_Matrix_View(P, &M, P->NbConstraints);
	return constraints_opt(&M, obj, denom, dir, opt, options);
    }
}

enum lp_result polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    if (options->lp_solver == BV_LP_POLYLIB)
	return PL_polyhedron_range(D, obj, denom, min, max, options);
    else if (options->lp_solver == BV_LP_GLPK)
	return glpk_polyhedron_range(D, obj, denom, min, max, options);
    else if (options->lp_solver == BV_LP_CDD)
	return cdd_polyhedron_range(D, obj, denom, min, max, options);
    else if (options->lp_solver == BV_LP_CDDF)
	return cddf_polyhedron_range(D, obj, denom, min, max, options);
    else if (options->lp_solver == BV_LP_PIP)
	return pip_polyhedron_range(D, obj, denom, min, max, options);
    else
	assert(0);
}
