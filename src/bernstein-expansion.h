/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include "polynomial.h"

#define DIGITS 10

string int2String(int n);
ex replaceVariablesInPolynomial(ex &poly, const GiNaC::exvector& V, ex &variables);
ex getBasis(unsigned int nbVert, matrix &A);
matrix getAiMatrix(unsigned int nbVert);
int bernsteinExpansion(matrix &P, ex &poly, const GiNaC::exvector& V,
		       unsigned int nbVert, unsigned int maxDegree,
		       const GiNaC::exvector& params);
bool linearCoefficients(lst coeffs, const GiNaC::exvector &Params);
void generateMaxConstraints(lst coeffs, const GiNaC::exvector& Params);
ex constantTerm(ex poly, const GiNaC::exvector &Params);
ex calculateLCM(ex poly, const GiNaC::exvector &Params, ex &cte);
extern "C" unsigned int checkConstraint(long long *M, unsigned int rows, unsigned int columns);
bool generateMaxConstraints(lst coeffs, const GiNaC::exvector &Params, unsigned int max);
bool generateMinConstraints(lst coeffs, const GiNaC::exvector &Params, unsigned int min);
bool generatePositiveNegativeConstraints(bool positive);
ex getMaxMinCoefficient1Param(lst coeffs, unsigned int maxDegree, ex Param, bool max, bool positive);
bool calculateDirection(bool max, bool positive, bool even);
void ex2longlongRow(long long *M, unsigned int row, polynomial &difference,
		    const GiNaC::exvector& Params);

ex powerMonomials(polynomial &poly, matrix &A, unsigned int nbVert
		  , unsigned int maxDegree, ex &basis);
lst getCoefficients(ex &maxDegreePolynomial, polynomial &expandedBasis
		    , unsigned int nbVert, matrix &A);
