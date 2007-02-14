#if defined(__cplusplus)
extern "C" {
#endif

#include <barvinok/evalue.h>
#include "verify.h"

/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

int check_poly(Polyhedron *S,Polyhedron *C,Enumeration *en,
	       int nparam, int pos, Value *z, const struct verify_options *options);
#if defined(__cplusplus)
}
#endif
