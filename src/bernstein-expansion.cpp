/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include <iostream>
#include <cln/cln.h>
#include <string>

#include "bernstein-expansion.h"

using namespace std;
using namespace GiNaC;

#define DIGITS 10

static std::string int2String(int n);
static unsigned int findMaxDegree(ex polynomial, const exvector& Vars);
static matrix getAiMatrix(unsigned int nbVert);
static ex getBasis(unsigned int nbVert, matrix &A);
static ex replaceVariablesInPolynomial(ex &poly, const exvector& V, ex &variables);
static ex powerMonomials(polynomial &poly, matrix &A, unsigned int nbVert,
			 unsigned int maxDegree, ex &basis);
static lst getCoefficients(ex &maxDegreePolynomial, polynomial &expandedBasis,
			   unsigned int nbVert, matrix &A);

exvector constructParameterVector(char **param_names, unsigned nbParams)
{
	exvector P(nbParams);
	for (int i = 0; i < nbParams; ++i) {
		P[i] = symbol(param_names[i]);
#ifdef DEBUG
		cout << "P: " << P[i] << endl;
#endif
	}
	return P;
}


exvector constructVariableVector(unsigned nbVariables, const char *prefix)
{
	exvector V(nbVariables);
	for (int i = 0; i < nbVariables; ++i) {
		V[i] = symbol(prefix + int2String(i));
#ifdef DEBUG
		cout << "V: " << V[i] << endl;
#endif
	}
	return V;
}


static numeric value2numeric(Value v)
{
    int sa = v[0]._mp_size;
    if (!sa)
	return 0;
    int abs_sa = sa < 0 ? -sa : sa;
    cln::cl_I res = 0;
    for (int i = abs_sa-1; i >= 0; --i) {
	res = res << GMP_LIMB_BITS;
	res = res + v[0]._mp_d[i];
    }
    return numeric(sa < 0 ? -res : res);
}


matrix domainVertices(Param_Polyhedron *PP, Param_Domain *Q, const exvector& params)
{
	Param_Vertices *V;
	unsigned nbVertices = 0;
	unsigned nbRows, nbColumns;
	int v;

	assert(PP->nbV > 0);
	nbRows = PP->V->Vertex->NbRows;
	nbColumns = PP->V->Vertex->NbColumns;

	FORALL_PVertex_in_ParamPolyhedron(V, Q, PP)
		++nbVertices;
	END_FORALL_PVertex_in_ParamPolyhedron;

	matrix VM(nbVertices, nbRows);		// vertices matrix

	v = 0;
	FORALL_PVertex_in_ParamPolyhedron(V, Q, PP)
		for (unsigned i = 0; i < nbRows; i++) {
			ex t;
			for (unsigned j = 0; j < nbColumns-2; j++)
				t += value2numeric(V->Vertex->p[i][j]) * params[j];
			t += value2numeric(V->Vertex->p[i][nbColumns-2]);
			t /= value2numeric(V->Vertex->p[i][nbColumns-1]);
#ifdef DEBUG
			cout << "T: " << t << endl;
#endif
			VM(v, i) = t;
		}
		++v;
	END_FORALL_PVertex_in_ParamPolyhedron;

	return VM;
}

/*
 * Do the Bernstein Expansion 
 *
 *	vert: vertices matrix
 *	poly: polynomial expression
 *	vars: vector of variables
 *	Params: vector of parameter
 */
lst bernsteinExpansion(const matrix &vert, ex poly, const exvector& vars,
		       const exvector& Params)
{
	unsigned maxDegree = findMaxDegree(poly, vars);
	matrix A = getAiMatrix(vert.rows());

#ifdef DEBUG
	cout << endl << "Polynomial: " << poly << endl << endl;
	cout << vert << endl << endl;
	cout << A << endl << endl;
#endif

	// obtain the variables value on the basis and replace
	ex substitutions = evalm(A * vert);
	ex polynom = replaceVariablesInPolynomial(poly, vars, substitutions);

#ifdef DEBUG
	cout << variables << endl << endl;
	cout << "Preliminar Expansion: " << polynom << endl << endl;
#endif

	ex basis = getBasis(vert.rows(), A);
	ex basisPowered = pow(basis, maxDegree);
	polynomial expandedBasis = basisPowered.expand();

#ifdef DEBUG
	cout << "Basis: " << basis << endl<< endl;
	cout << "Powered Basis: " << basisPowered << endl;
	cout << "Expanded Basis: " << expandedBasis << endl << endl;
#endif

	// monomials to n degree
	polynomial p(polynom);
	ex maxDegreePolynomial = powerMonomials(p, A, vert.rows(), maxDegree, basis);

	// get the coefficients
	return getCoefficients(maxDegreePolynomial, expandedBasis, vert.rows(), A);
}

/*
 * Construct A_i matrix
 *
 *	nbVert: number of vertices
 */
matrix getAiMatrix(unsigned int nbVert)
{
	matrix A(1, nbVert); 		// a_i matrix
	for(unsigned int i = 0; i < nbVert; i++) {
		A(0,i) = symbol("a" + int2String(i));
#ifdef DEBUG
		cout << "A: " << A(0,i) << endl;
#endif
	}
	return A;
}


/*
 * Construct the basis
 *
 *	A: a_i matrix
 *	nbVert: number of vertices
 */
ex getBasis(unsigned int nbVert, matrix &A)
{
	ex basis;
	for(unsigned int i = 0; i < nbVert; i++) {
		basis += A(0,i);
	}
	return basis;
}


/*
 * Finds the less than maxDegree monomials and multiply them
 * by the basis
 *
 *	polynomial: original polynomial
 *	A: a_i matrix
 *	nbVert: number of vertices
 *	maxDegree: maximal polynomial multidegree
 *	basis: basis of the polytope
 */
ex powerMonomials(polynomial &poly, matrix &A, unsigned int nbVert
		  , unsigned int maxDegree, ex &basis)
{
	ex maxDegreePolynomial;
#ifdef DEBUG
	cout << "- Degree --------------------------------------" << endl;
#endif
	for (size_t i = 0; i != poly.nbTerms(); ++i) {
		unsigned int degree = 0;		// degree of the monomial

		for(unsigned int j = 0; j < nbVert; j++) {
			degree += poly.term(i).degree(A(0,j));
		}
#ifdef DEBUG
		cout << poly.term(i) << " Degree: " << degree;
#endif
		if(degree < maxDegree) {
			ex degreeUp = poly.term(i) * pow(basis, maxDegree - degree);
#ifdef DEBUG
			cout << "   --> New Term: " <<  degreeUp.expand();
#endif
			maxDegreePolynomial += degreeUp.expand();
		} else {
			maxDegreePolynomial += poly.term(i);
		}
#ifdef DEBUG
		cout << endl << "-----------------------------------------------" << endl;
#endif

	}
#ifdef DEBUG
	cout << endl << "Final Expansion: " << maxDegreePolynomial << endl << endl;
#endif
	return maxDegreePolynomial;
}



/*
 * Finds and prints the coefficients of the polynomial
 *
 *	maxDegreePolynomial: polynomial with monomials of degree n
 *	expandedBasis: basis powered n
 *	nbVert: number of vertices
 *	A: a_i matrix
 */
lst getCoefficients(ex &maxDegreePolynomial, polynomial &expandedBasis
		       , unsigned int nbVert, matrix &A)
{
	lst coefficients;

	for (size_t i = 0; i != expandedBasis.nbTerms(); ++i) {
		ex coeff = maxDegreePolynomial;
		ex vars = 1;
		for(unsigned int j = 0; j < nbVert; j++) {
			unsigned int deg = expandedBasis.term(i).degree(A(0,j));
			if(deg > 0) {
				coeff = coeff.coeff(A(0,j), deg);
				vars *= pow(A(0,j), deg);
			}
		}
		coefficients.append(coeff / (expandedBasis.term(i)/vars));
	}
	return coefficients.sort().unique();
}


// replace the variables in the polynomial
ex replaceVariablesInPolynomial(ex &poly, const exvector &vars, ex &substitutions)
{
	lst replace;

	for(unsigned int i = 0; i < vars.size(); i++) {
#ifdef DEBUG
		cout << "Replacing: " << vars[i] << " by " << substitutions[i] << endl;
#endif
		replace.append(vars[i] == substitutions[i]);
	}
	ex polyRepl = poly.subs(replace);

	return(polyRepl.expand());
}


/* Converts int n to string */
string int2String(int n)
{
	char numeroVariable[DIGITS];
	snprintf(numeroVariable, DIGITS, "%d", n);
	string nroV(numeroVariable);

	return nroV;
}


/* 
 * Find the maximum multi-degree of the polinomial
 *
 *	polynomial: polynomial
 *	Vars: variables matrix
 *	nbVar: number of variables
 */
unsigned int findMaxDegree(ex polynomial, const exvector& Vars)
{
	unsigned int max = polynomial.degree(Vars[0]);

	for(unsigned int i = 0; i < Vars.size(); i++) {
		unsigned int degree = polynomial.degree(Vars[i]);
		if(max < degree) {
			max = degree;
		}
	}
	return max;
}
