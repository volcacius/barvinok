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
#include "bernstein.h"
#include "bernstein-expansion.h"


/* 
 * Convert a matrix from long long C-style matrix to Ginac expression
 *
 *	m: polynomial matrix
 *	nbRows, nbColumns: number of rows, columns
 *	Vars: variables of the polynomial
 *	nbVar: number of variables
 */
ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns, 
		     const exvector& Vars)
{
	ex p;

	for(unsigned int i = 0; i < nbRows; i++) {
		ex t;
		for(unsigned int j = 0; j < Vars.size(); j++) {
			// TODO: loosing precision to long int
			long int val = (long int) m[i*nbColumns+j];
			if(val != 0) {
				t += pow(Vars[j], (long int) m[i*nbColumns+j]);
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
ex polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns,
			  long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns,
			  const exvector& vars, const exvector& params)
{
	unsigned int nbVariables = nbColumns;
	unsigned int nbCoefficients = nbRows;
	ex p;

	for(unsigned int i = 0; i < nbCoefficients; i++) {
#ifdef DEBUG
		cout << "-----------------------" << endl;
#endif
		ex c;
		ex t;

		c = convertPolynomial(llPolynomialCoefficients[i], llRows[i], llColumns[i], params);
#ifdef DEBUG
		cout << "Coeff[i]: " << c << endl;
#endif

		for(unsigned int j = 0; j < nbVariables; j++) {
			// TODO: loosing precision to long int
			long int val = (long int) m[i*nbColumns+j];
			if(val != 0) {
				t += pow(vars[j], (long int) m[i*nbColumns+j]);
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

#ifdef DEBUG
	cout << "Polynomial: " << p << endl;
#endif

	return p;
}
