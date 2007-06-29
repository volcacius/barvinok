#include <barvinok/options.h>
#include <barvinok/util.h>
#include "polysign.h"
#include "config.h"

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
#endif

enum order_sign polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    if (options->lp_solver == BV_LP_POLYLIB)
	return PL_polyhedron_affine_sign(D, T, options);
    else if (options->lp_solver == BV_LP_CDD)
	return cdd_polyhedron_affine_sign(D, T, options);
    else if (options->lp_solver == BV_LP_CDDF)
	return cddf_polyhedron_affine_sign(D, T, options);
    else
	assert(0);
}

enum lp_result polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    if (options->lp_solver == BV_LP_POLYLIB)
	return PL_polyhedron_range(D, obj, denom, min, max, options);
    else if (options->lp_solver == BV_LP_CDD)
	return cdd_polyhedron_range(D, obj, denom, min, max, options);
    else if (options->lp_solver == BV_LP_CDDF)
	return cddf_polyhedron_range(D, obj, denom, min, max, options);
    else
	assert(0);
}
