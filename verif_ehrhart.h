#if defined(__cplusplus)
extern "C" {
#endif

#include <barvinok/evalue.h>
#include "verify.h"

/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

Polyhedron *check_poly_context_scan(Polyhedron *C, struct verify_options *options);
void check_poly_init(Polyhedron *C, struct verify_options *options);
int check_poly(Polyhedron *S, Polyhedron *CS, evalue *EP, int exist,
	       int nparam, int pos, Value *z, const struct verify_options *options);
#if defined(__cplusplus)
}
#endif
