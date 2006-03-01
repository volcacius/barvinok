/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"
#include "polynomial.h"

#define DIGITS 10

std::string int2String(int n);
GiNaC::matrix domainVertices(Param_Polyhedron *PP, Param_Domain *Q, 
			     const GiNaC::exvector& params);
GiNaC::lst bernsteinExpansion(const GiNaC::matrix &vert, GiNaC::ex poly, 
			      const GiNaC::exvector& vars,
			      const GiNaC::exvector& params);
int getMaxMinCoefficient(Polyhedron *VD, GiNaC::lst coeffs, 
			const GiNaC::exvector& Params);
void generateMaxConstraints(GiNaC::lst coeffs, const GiNaC::exvector& Params);
GiNaC::ex constantTerm(GiNaC::ex poly, const GiNaC::exvector &Params);
GiNaC::ex calculateLCM(GiNaC::ex poly, const GiNaC::exvector &Params, GiNaC::ex &cte);
bool generateMaxConstraints(Polyhedron *VD, GiNaC::lst coeffs, 
			    const GiNaC::exvector &Params, unsigned int max);
bool generateMinConstraints(Polyhedron *VD, GiNaC::lst coeffs, 
			    const GiNaC::exvector &Params, unsigned int min);
bool generatePositiveNegativeConstraints(Polyhedron *VD, bool positive);
GiNaC::ex getMaxMinCoefficient1Param(GiNaC::lst coeffs, unsigned int maxDegree, 
				     GiNaC::ex Param, bool max, bool positive);
bool calculateDirection(bool max, bool positive, bool even);
void ex2longlongRow(long long *M, unsigned int row, polynomial &difference,
		    const GiNaC::exvector& Params);

