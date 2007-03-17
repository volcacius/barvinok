/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include <ginac/ginac.h>
#include <gmp.h>
extern "C" {
#define matrix polylib_matrix
#define polynomial polylib_polynomial
#include <polylib/polylibgmp.h>
#undef matrix
#undef polynomial
#undef value_compare
#undef divide
}

namespace bernstein {

GiNaC::numeric value2numeric(const Value v);
void numeric2value(const GiNaC::numeric& n, Value& v);
GiNaC::exvector constructParameterVector(const char * const *param_names, 
					 unsigned nbParams);
GiNaC::exvector constructVariableVector(unsigned nbVariables, const char *prefix);
GiNaC::ex replaceVariablesInPolynomial(const GiNaC::ex &poly,
				       const GiNaC::exvector& V,
				       const GiNaC::ex &variables);
GiNaC::matrix domainVertices(Param_Polyhedron *PP, Param_Domain *Q, 
			     const GiNaC::exvector& params);
GiNaC::lst bernsteinExpansion(const GiNaC::matrix& vert, const GiNaC::ex& poly, 
			      const GiNaC::exvector& vars,
			      const GiNaC::exvector& params);

}
