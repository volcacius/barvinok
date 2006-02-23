/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#define MAXRAYS 1000

#if defined(__cplusplus)
extern "C" {
#endif

extern int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns
				     , long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns
				     , unsigned int nbParams, char **param_values);
extern int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbParams, char **param_values);

long long *matrix2longlong(Matrix *M);

unsigned  checkConstraint(long long *M, unsigned int rows, unsigned int columns);
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);
void Matrix_Set(int r, int c, int v, Matrix * A);

void printDomain(Param_Domain *Q, Polyhedron *B, char **param_name);

#if defined(__cplusplus)
}
#endif
