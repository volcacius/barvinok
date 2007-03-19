#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

void Param_Polyhedron_Scale_Integer_Slow(Param_Polyhedron *PP, Polyhedron **P,
					 Value *det, unsigned MaxRays);

#if defined(__cplusplus)
}
#endif
