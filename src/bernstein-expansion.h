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
}
#include "polynomial.h"

#define DIGITS 10

string int2String(int n);
ex replaceVariablesInPolynomial(ex &poly, const GiNaC::exvector& V, ex &variables);
ex getBasis(unsigned int nbVert, matrix &A);
matrix getAiMatrix(unsigned int nbVert);
bool linearCoefficients(lst coeffs, const GiNaC::exvector &Params);
void generateMaxConstraints(lst coeffs, const GiNaC::exvector& Params);
ex constantTerm(ex poly, const GiNaC::exvector &Params);
ex calculateLCM(ex poly, const GiNaC::exvector &Params, ex &cte);
bool generateMaxConstraints(Polyhedron *VD, lst coeffs, 
GiNaC::lst bernsteinExpansion(GiNaC::matrix &P, GiNaC::ex &poly, 
			      const GiNaC::exvector& V,
			      unsigned int nbVert, const GiNaC::exvector& params);
int getMaxMinCoefficient(Polyhedron *VD, GiNaC::lst coeffs, 
			const GiNaC::exvector& Params);
			    const GiNaC::exvector &Params, unsigned int max);
bool generateMinConstraints(Polyhedron *VD, lst coeffs, 
			    const GiNaC::exvector &Params, unsigned int min);
bool generatePositiveNegativeConstraints(Polyhedron *VD, bool positive);
ex getMaxMinCoefficient1Param(lst coeffs, unsigned int maxDegree, ex Param, bool max, bool positive);
bool calculateDirection(bool max, bool positive, bool even);
void ex2longlongRow(long long *M, unsigned int row, polynomial &difference,
		    const GiNaC::exvector& Params);

ex powerMonomials(polynomial &poly, matrix &A, unsigned int nbVert
		  , unsigned int maxDegree, ex &basis);
lst getCoefficients(ex &maxDegreePolynomial, polynomial &expandedBasis
		    , unsigned int nbVert, matrix &A);
