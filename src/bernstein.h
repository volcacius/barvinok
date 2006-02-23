/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#define MAXRAYS 1000

extern int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns
				     , long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns
				     , unsigned int nbParams, char **param_values);
extern int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbParams, char **param_values);

Matrix *readPolynomial(unsigned int nbVariables, unsigned int nbParams, char **param_name);
long long *matrix2longlong(Matrix *M);

unsigned  checkConstraint(long long *M, unsigned int rows, unsigned int columns);
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);
void Matrix_Set(int r, int c, int v, Matrix * A);

