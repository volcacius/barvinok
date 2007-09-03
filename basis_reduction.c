#include <assert.h>
#include <barvinok/basis_reduction.h>
#include <barvinok/options.h>
#include "config.h"

#ifndef HAVE_LIBGLPK
Matrix *glpk_Polyhedron_Reduced_Basis(Polyhedron *P)
{
    assert(0);
}
#endif

#ifndef HAVE_LIBCDDGMP
Matrix *cdd_Polyhedron_Reduced_Basis(Polyhedron *P)
{
    assert(0);
}
#endif

Matrix *Polyhedron_Reduced_Basis(Polyhedron *P, struct barvinok_options *options)
{
    if (options->gbr_lp_solver == BV_GBR_GLPK)
	return glpk_Polyhedron_Reduced_Basis(P);
    else if (options->gbr_lp_solver == BV_GBR_CDD)
	return cdd_Polyhedron_Reduced_Basis(P);
    else
	assert(0);
}
