/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#define matrix polylib_matrix
#define polynomial polylib_polynomial
#include <polylib/polylibgmp.h>
#undef matrix
#undef polynomial
#undef value_compare
#undef divide

#define MAXRAYS 1000

unsigned  checkConstraint(Polyhedron *VD, long long *M, 
			  unsigned int rows, unsigned int columns);
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);

#if defined(__cplusplus)
}
#endif
