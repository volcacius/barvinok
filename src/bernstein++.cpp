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



/* 
 * Convert vertices from long long C-style matrix to Ginac-style matrix 
 *
 *	m: vertex matrix
 *	nbRows, nbColumns: number of rows, columns
 *	nbVertices: number of vertices in the matrix
 *	param_values: parameters names
 *	nbParams: number of parameters
 *
 */
extern "C" int verticesConvert(long long *m, unsigned int nbRows, unsigned int nbColumns
			       , unsigned int nbVertices, char **param_values, unsigned int nbParams)
{
	matrix V(nbVertices, nbRows);		// vertices matrix

#ifdef DEBUG
	cout << "polynomial: " << polynomial << endl;
	cout << "Vertices: " << nbVertices << endl;
	cout << "Col: " << nbColumns << endl;
	cout << "Row per vertice: " << nbRows << endl;
	cout << "Size: " << nbVertices * nbColumns * nbRows<<endl;
	cout << "P: " << P << endl;
#endif

	for(unsigned int v = 0; v < nbVertices; v++) {
		for(unsigned int i = 0; i < nbRows; i++) {
			ex t;
			for(unsigned int j = 0; j < nbColumns-2; j++) {
				// TODO: loosing precision to long int
				t += P(0, j) *  (long int) m[v*nbColumns*nbRows+i*nbColumns+j];
			}
			t += (long int) m[v*nbColumns*nbRows+i*nbColumns+nbColumns-2];
			t /= (long int) m[v*nbColumns*nbRows+i*nbColumns+nbColumns-1];
#ifdef DEBUG
			cout << "T: " << t << endl;
#endif
			V(v, i) = t;
		}
	}
#ifdef DEBUG
	cout << "Matrix V: " << V << endl;
#endif

	// do the expansion
	bernsteinExpansion(V, polynomial, Vars, nbVar, nbVertices, findMaxDegree(polynomial, Vars, nbVar), P, nbParams);

	return 0;
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
