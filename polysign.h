#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include "lexmin.h"

enum order_sign { order_lt, order_le, order_eq, order_ge, order_gt, order_unknown,
		  order_undefined };

/* Returns the sign of the affine function specified by T on the polyhedron D */
enum order_sign polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);
enum order_sign PL_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);
enum order_sign cdd_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);
enum order_sign cddf_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options);

enum lp_result { lp_ok = 0, lp_unbounded, lp_empty };

/* Computes min and max of obj/denom over D, with min rounded up
 * and max rounded down to an integer.  obj has D->Dimension+1
 * elements, with the constant term in position D->Dimension.
 */
enum lp_result polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options);
enum lp_result PL_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options);
enum lp_result glpk_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options);
enum lp_result cdd_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options);
enum lp_result cddf_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options);
#if defined(__cplusplus)
}
#endif
