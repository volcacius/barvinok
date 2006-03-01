/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#include <gmp.h>
#include "polylib++.h"

#if defined(__cplusplus)
extern "C" {
#endif

#define MAXRAYS 1000

unsigned  checkConstraint(Polyhedron *VD, long long *M, 
			  unsigned int rows, unsigned int columns);
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);

#if defined(__cplusplus)
}
#endif
