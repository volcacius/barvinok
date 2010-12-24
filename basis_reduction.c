#include <assert.h>
#include <isl_set_polylib.h>
#include <barvinok/basis_reduction.h>
#include <barvinok/options.h>
#include "config.h"

#ifndef HAVE_LIBGLPK
Matrix *glpk_Polyhedron_Reduced_Basis(Polyhedron *P,
				     struct barvinok_options *options)
{
    assert(0);
}
#endif

#ifndef HAVE_LIBCDDGMP
Matrix *cdd_Polyhedron_Reduced_Basis(Polyhedron *P,
				     struct barvinok_options *options)
{
    assert(0);
}
#endif

Matrix *isl_Polyhedron_Reduced_Basis(Polyhedron *P,
				     struct barvinok_options *options)
{
	int i, j;
	isl_int v;
	isl_ctx *ctx;
	isl_dim *dim;
	int nvar = P->Dimension;
	isl_basic_set *bset;
	isl_mat *basis;
	Matrix *M;
	int isl_gbr_only_first = options->isl->gbr_only_first;

	options->isl->gbr_only_first = options->gbr_only_first;
	ctx = isl_ctx_alloc_with_options(barvinok_options_arg, options);
	assert(ctx);

	dim = isl_dim_set_alloc(ctx, 0, nvar);
	bset = isl_basic_set_new_from_polylib(P, dim);

	basis = isl_basic_set_reduced_basis(bset);
	isl_basic_set_free(bset);

	M = Matrix_Alloc(nvar, nvar);

	isl_int_init(v);
	for (i = 0; i < nvar; ++i)
		for (j = 0; j < nvar; ++j) {
			isl_mat_get_element(basis, 1 + i, 1 + j, &v);
			isl_int_get_gmp(v, M->p[i][j]);
		}
	isl_int_clear(v);

	isl_mat_free(basis);

	isl_ctx_free(ctx);

	options->isl->gbr_only_first = isl_gbr_only_first;

	return M;
}

Matrix *Polyhedron_Reduced_Basis(Polyhedron *P, struct barvinok_options *options)
{
    if (options->gbr_lp_solver == BV_GBR_GLPK)
	return glpk_Polyhedron_Reduced_Basis(P, options);
    else if (options->gbr_lp_solver == BV_GBR_CDD)
	return cdd_Polyhedron_Reduced_Basis(P, options);
    else if (options->gbr_lp_solver == BV_GBR_ISL)
	return isl_Polyhedron_Reduced_Basis(P, options);
    else
	assert(0);
}
