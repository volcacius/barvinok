/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include "polynomial.h"

#define DIGITS 10

string int2String(int n);
ex replaceVariablesInPolynomial(ex &poly, matrix &V, unsigned int nbVar, ex &variables);
ex getBasis(unsigned int nbVert, matrix &A);
matrix getAiMatrix(unsigned int nbVert);
int bernsteinExpansion(matrix &P, ex &poly, matrix &V, unsigned int nbVar
		       , unsigned int nbVert, unsigned int maxDegree
		       , matrix &Params, unsigned int nbParams);
bool linearCoefficients(lst coeffs, matrix &Params, unsigned int nbParams);
void generateMaxConstraints(lst coeffs, matrix &Params, unsigned int nbParams);
ex constantTerm(ex poly, matrix &Params, unsigned int nbParams);
ex calculateLCM(ex poly, matrix &Params, unsigned int nbParams, ex &cte);
extern "C" unsigned int checkConstraint(long long *M, unsigned int rows, unsigned int columns);
bool generateMaxConstraints(lst coeffs, matrix &Params, unsigned int nbParams, unsigned int max);
bool generateMinConstraints(lst coeffs, matrix &Params, unsigned int nbParams, unsigned int min);
bool generatePositiveNegativeConstraints(matrix &Params, bool positive);
ex getMaxMinCoefficient1Param(lst coeffs, unsigned int maxDegree, matrix &Params, bool max, bool positive);
bool calculateDirection(bool max, bool positive, bool even);
void ex2longlongRow(long long *M, unsigned int row, polynomial &difference, matrix &Params, unsigned int nbParams);

ex powerMonomials(polynomial &poly, matrix &A, unsigned int nbVert
		  , unsigned int maxDegree, ex &basis);
lst getCoefficients(ex &maxDegreePolynomial, polynomial &expandedBasis
		    , unsigned int nbVert, matrix &A);
