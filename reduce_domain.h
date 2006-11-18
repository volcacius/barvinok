#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

Polyhedron *reduce_domain(Polyhedron *D, Matrix *CT, Polyhedron *CEq,
			  Polyhedron **fVD, int nd, unsigned MaxRays);

#if defined(__cplusplus)
}
#endif
