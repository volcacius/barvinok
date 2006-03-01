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

GiNaC::numeric value2numeric(Value v);
GiNaC::exvector constructParameterVector(char **param_names, unsigned nbParams);
GiNaC::exvector constructVariableVector(unsigned nbVariables, const char *prefix);
GiNaC::matrix domainVertices(Param_Polyhedron *PP, Param_Domain *Q, 
			     const GiNaC::exvector& params);
GiNaC::lst bernsteinExpansion(const GiNaC::matrix &vert, GiNaC::ex poly, 
			      const GiNaC::exvector& vars,
			      const GiNaC::exvector& params);

}
