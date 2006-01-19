/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#define MAXRAYS 1000

extern int verticesConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbVertices
			   , char **param_values, unsigned int nbParams);

extern int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns
				     , long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns
				     , unsigned int nbParams, char **param_values);
extern int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbParams, char **param_values);

Matrix *readPolynomial(unsigned int nbVariables, unsigned int nbParams, char **param_name);
long long *matrix2longlong(Matrix *M);

unsigned  checkConstraint(long long *M, unsigned int rows, unsigned int columns);
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);
void Matrix_Set(int r, int c, int v, Matrix * A);
void convertMatrix(long long *matrix, Param_Polyhedron *PP, Param_Domain *Q, char **param_name
		   , unsigned int nbColumns, unsigned int nbRows);

