/* 
 *	Bernstein Expansion
 *
 *	- c to c++ functions
 */

#include <string.h>
#include <ginac/ginac.h>
#include <gmp.h>
extern "C" {
#define matrix polylib_matrix
#define polynomial polylib_polynomial
#include <polylib/polylibgmp.h>
#undef matrix
#undef polynomial
}

unsigned int findMaxDegree(GiNaC::ex &polynomial, GiNaC::exvector Vars);

matrix getParameterMatrix(unsigned int nbParams, char **param_values);
GiNaC::exvector constructParameterVector(char **param_names, unsigned nbParams);
void getVariablesMatrix(unsigned int nbVariables);
GiNaC::exvector constructVariableVector(unsigned nbVariables, const char *prefix);
std::string int2String(int n);
GiNaC::ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns,
			    const GiNaC::exvector& params);
Matrix *readPolynomial(unsigned int nbVariables, const GiNaC::exvector& params);

int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns, 
			  long long **llPolynomialCoefficients, unsigned int *llRows, 
			  unsigned int *llColumns, const GiNaC::exvector& params);
int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns);

extern "C" {

void doExpansion(Param_Polyhedron *PP, Param_Domain *Q, GiNaC::exvector params);

}
