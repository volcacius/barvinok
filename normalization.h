#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

Matrix *standard_constraints(Polyhedron *P, unsigned nparam, int *rows_p,
			     Matrix **T);

#if defined(__cplusplus)
}
#endif
