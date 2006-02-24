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

extern int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns
				     , long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns
				     , unsigned int nbParams, char **param_values);
extern int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbParams, char **param_values);

long long *matrix2longlong(Matrix *M);

unsigned  checkConstraint(Polyhedron *VD, long long *M, 
			  unsigned int rows, unsigned int columns);
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);
void Matrix_Set(int r, int c, int v, Matrix * A);

void printDomain(Param_Domain *Q, Polyhedron *B, char **param_name);

#if defined(__cplusplus)
}
#endif
