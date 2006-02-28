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

GiNaC::exvector constructParameterVector(char **param_names, unsigned nbParams);
GiNaC::exvector constructVariableVector(unsigned nbVariables, const char *prefix);
std::string int2String(int n);
GiNaC::ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns,
			    const GiNaC::exvector& params);

GiNaC::ex polyConvertParameters(long long *m, unsigned int nbRows, 
			  unsigned int nbColumns, 
			  long long **llPolynomialCoefficients, unsigned int *llRows, 
			  unsigned int *llColumns, const GiNaC::exvector& vars,
			  const GiNaC::exvector& params);
