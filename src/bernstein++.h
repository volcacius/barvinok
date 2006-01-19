/* 
 *	Bernstein Expansion
 *
 *	- c to c++ functions
 */

unsigned int findMaxDegree(ex &polynomial, matrix &Vars, unsigned int nbVar);

int bernsteinExpansion(matrix &P, ex &poly, matrix &V, unsigned int nbVar
		       , unsigned int nbVert, unsigned int maxDegree
		       , matrix &Params, unsigned int nbParams);
matrix getParameterMatrix(unsigned int nbParams, char **param_values);
void getVariablesMatrix(unsigned int nbVariables);
string int2String(int n);
ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns, matrix &Vars, unsigned int nbVariables);

extern "C" int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns
				     , long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns
				     , unsigned int nbParams, char **param_values);
extern "C" int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbParams, char **param_values);
