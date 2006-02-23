/* 
 *	Bernstein Expansion
 *
 *	- c to c++ functions
 */



#include <iostream>
#include <string>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

#include "bernstein++.h"


/* Global variables to send to bernstein_expansion */
ex polynomial;		// polynomial

matrix Vars;		// variables matrix
unsigned int nbVar;		// number of variables

matrix P;			// parameters matrix


/*
 * Construct Parameter matrix
 *
 *	nbParams: number of parameters
 *	params_values: names of the parameters
 */
matrix getParameterMatrix(unsigned int nbParams, char **param_values)
{
	matrix P(1, nbParams); 		// parameters matrix
	symbol *pSym;

	for(unsigned int i = 0; i < nbParams; i++) {
		pSym = new symbol(param_values[i]);
#ifdef DEBUG
		cout << "P: " << *pSym << endl;
#endif
		P(0,i) = *pSym;
		delete pSym;
	}
	return P;
}



/* Given a domain, converts the vertices matrix to an long long matrix and then
   sends it to ginac functions */
void doExpansion(Param_Polyhedron *PP, Param_Domain *Q, unsigned int nb_param,
		 char **param_name)
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
			for (unsigned j = 0; j < nbColumns-2; j++) {
				// TODO: losing precision to long int
				t += VALUE_TO_LONG(V->Vertex->p[i][j]) * P(0, j);
			}
			t += VALUE_TO_LONG(V->Vertex->p[i][nbColumns-2]);
			t /= VALUE_TO_LONG(V->Vertex->p[i][nbColumns-1]);
#ifdef DEBUG
			cout << "T: " << t << endl;
#endif
			VM(v, i) = t;
		}
		++v;
	END_FORALL_PVertex_in_ParamPolyhedron;

	bernsteinExpansion(VM, polynomial, Vars, nbVar, nbVertices, findMaxDegree(polynomial, Vars, nbVar), P, nb_param);
}


/* 
 * Find the maximum multi-degree of the polinomial
 *
 *	polynomial: polynomial
 *	Vars: variables matrix
 *	nbVar: number of variables
 */
unsigned int findMaxDegree(ex &polynomial, matrix &Vars, unsigned int nbVar)
{
	unsigned int max = polynomial.degree(Vars(0, 0));

	for(unsigned int i = 0; i < nbVar; i++) {
		unsigned int degree = polynomial.degree(Vars(0,i));
		if(max < degree) {
			max = degree;
		}
	}
	return max;
}


/*
 * Construct Variables matrix
 *
 *	nbVariables: number of variables
 */
void getVariablesMatrix(unsigned int nbVariables)
{
	symbol *pSym;

	for(unsigned int i = 0; i < nbVariables; i++) {
		pSym = new symbol("v" + int2String(i));
#ifdef DEBUG
		cout << "V: " << *pSym << endl;
#endif
		Vars(0,i) = *pSym;
		delete pSym;
	}
}


/* 
 * Convert polynomial from long long C-style matrix to Ginac expression
 *
 *	m: polynomial matrix
 *	nbRows, nbColumns: number of rows, columns
 */
extern "C" int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns, unsigned int nbParams, char **param_values)
{
	unsigned int nbVariables = nbColumns - 2;
	ex p;

	P = getParameterMatrix(nbParams, param_values);	// parameters matrix

	// setting global variables
	Vars = matrix(1, nbVariables); 		// variables matrix
	nbVar = nbVariables;

	getVariablesMatrix(nbVariables);
       
	p = convertPolynomial(m, nbRows, nbColumns, Vars, nbVar);
	polynomial = p;
	return 0;
}


/* 
 * Convert a matrix from long long C-style matrix to Ginac expression
 *
 *	m: polynomial matrix
 *	nbRows, nbColumns: number of rows, columns
 *	Vars: variables of the polynomial
 *	nbVar: number of variables
 */
ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns, matrix &Vars, unsigned int nbVariables)
{
	ex p;

	for(unsigned int i = 0; i < nbRows; i++) {
		ex t;
		for(unsigned int j = 0; j < nbVariables; j++) {
			// TODO: loosing precision to long int
			long int val = (long int) m[i*nbColumns+j];
			if(val != 0) {
				t += pow(Vars(0, j), (long int) m[i*nbColumns+j]);
#ifdef DEBUG
				cout << "T: " << t << endl;
#endif
			}
		}
		if(t == 0) {
			t += (long int) m[i*nbColumns+nbColumns-2];
		} else {
			t *= (long int) m[i*nbColumns+nbColumns-2];
		}
#ifdef DEBUG
		cout << "T: " << t << endl;
#endif
		t /= (long int) m[i*nbColumns+nbColumns-1];
#ifdef DEBUG
		cout << "T: " << t << endl;
		cout << endl;
#endif
		p += t;
	}
	return p;
}


/* 
 * Convert polynomial and coefficients from long long C-style matrix to Ginac expression
 *
 *	m: polynomial matrix
 *	nbRows, nbColumns: number of rows, columns
 *	llPolynomialCoefficients: coefficients matrices
 *	llRows: coefficients number of rows
 *	llColumns: coefficients number of columns
 */
extern "C" int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns
				     , long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns
				     , unsigned int nbParams, char **param_values)
{
	unsigned int nbVariables = nbColumns;
	unsigned int nbCoefficients = nbRows;
	ex p;

	P = getParameterMatrix(nbParams, param_values);	// parameters matrix

	// setting global variables
	Vars = matrix(1, nbVariables); 		// variables matrix
	nbVar = nbVariables;
	getVariablesMatrix(nbVariables);
       
	for(unsigned int i = 0; i < nbCoefficients; i++) {
#ifdef DEBUG
		cout << "-----------------------" << endl;
#endif
		ex c;
		ex t;

		c = convertPolynomial(llPolynomialCoefficients[i], llRows[i], llColumns[i], P, nbParams);
#ifdef DEBUG
		cout << "Coeff[i]: " << c << endl;
#endif

		for(unsigned int j = 0; j < nbVariables; j++) {
			// TODO: loosing precision to long int
			long int val = (long int) m[i*nbColumns+j];
			if(val != 0) {
				t += pow(Vars(0, j), (long int) m[i*nbColumns+j]);
#ifdef DEBUG
				cout << "T: " << t << endl;
#endif
			}
		}
#ifdef DEBUG
		cout << "T: " << t << endl;
#endif
		p += (t * c);
	}
	polynomial = p;

#ifdef DEBUG
	cout << "Polynomial: " << polynomial << endl;
#endif

	return 0;
}
