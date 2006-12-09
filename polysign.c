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
#endif

enum order_sign polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    if (options->lexmin_polysign == BV_LEXMIN_POLYSIGN_POLYLIB)
	return PL_polyhedron_affine_sign(D, T, options);
    else if (options->lexmin_polysign == BV_LEXMIN_POLYSIGN_CDD)
	return cdd_polyhedron_affine_sign(D, T, options);
    else if (options->lexmin_polysign == BV_LEXMIN_POLYSIGN_CDDF)
	return cddf_polyhedron_affine_sign(D, T, options);
    else
	assert(0);
}
