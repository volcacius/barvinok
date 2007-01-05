#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include "lexmin.h"

enum order_sign { order_lt, order_le, order_eq, order_ge, order_gt, order_unknown,
		  order_undefined };

/* Returns the sign of the affine function specified by T on the polyhedron D */
enum order_sign polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct lexmin_options *options);
enum order_sign polylib_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);
enum order_sign cdd_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);
enum order_sign cddf_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);

#if defined(__cplusplus)
}
#endif
