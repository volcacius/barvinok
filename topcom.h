#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

Matrix *standard_constraints(Polyhedron *P, unsigned nparam, int *rows_p,
			     Matrix **T);
Param_Polyhedron *TC_P2PP(Polyhedron *P, Polyhedron *C,
			  struct barvinok_options *options);

#if defined(__cplusplus)
}
#endif
