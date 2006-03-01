/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"
#include "polynomial.h"

GiNaC::exvector constructParameterVector(char **param_names, unsigned nbParams);
GiNaC::exvector constructVariableVector(unsigned nbVariables, const char *prefix);
GiNaC::matrix domainVertices(Param_Polyhedron *PP, Param_Domain *Q, 
			     const GiNaC::exvector& params);
GiNaC::lst bernsteinExpansion(const GiNaC::matrix &vert, GiNaC::ex poly, 
			      const GiNaC::exvector& vars,
			      const GiNaC::exvector& params);
