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

unsigned int findMaxDegree(GiNaC::ex &polynomial, GiNaC::matrix &Vars, 
			   unsigned int nbVar);

int bernsteinExpansion(GiNaC::matrix &P, GiNaC::ex &poly, GiNaC::matrix &V, 
		       unsigned int nbVar, unsigned int nbVert, 
		       unsigned int maxDegree, GiNaC::matrix &Params, 
		       unsigned int nbParams);
matrix getParameterMatrix(unsigned int nbParams, char **param_values);
void getVariablesMatrix(unsigned int nbVariables);
std::string int2String(int n);
GiNaC::ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns, matrix &Vars, unsigned int nbVariables);

extern "C" {

int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns, 
			  long long **llPolynomialCoefficients, unsigned int *llRows, 
			  unsigned int *llColumns, unsigned int nbParams, 
			  char **param_values);
int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, 
	        unsigned int nbParams, char **param_values);
void doExpansion(Param_Polyhedron *PP, Param_Domain *Q, unsigned int nb_param, 
		 char **param_name);

}
