#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

Polyhedron *reduce_domain(Polyhedron *D, Matrix *CT, Polyhedron *CEq,
			  Polyhedron **fVD, int nd, unsigned MaxRays);

#if defined(__cplusplus)
}
#endif
