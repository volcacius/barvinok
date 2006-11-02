#if defined(__cplusplus)
extern "C" {
#endif

#include <barvinok/evalue.h>

/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

extern Value Min, Max;

extern char **params;

extern int st;

int check_poly(Polyhedron *S,Polyhedron *C,Enumeration *en,
	       int nparam, int pos, Value *z, int print_all);
#if defined(__cplusplus)
}
#endif
