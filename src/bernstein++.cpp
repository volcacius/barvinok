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


/* Global variables to send to bernstein_expansion */
ex polynomial;		// polynomial

exvector Vars;		// variables matrix


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



/* Given a domain, converts the vertices matrix to an long long matrix and then
   sends it to ginac functions */
void doExpansion(Param_Polyhedron *PP, Param_Domain *Q, exvector params)
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
				t += VALUE_TO_LONG(V->Vertex->p[i][j]) * params[j];
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

	bernsteinExpansion(VM, polynomial, Vars, nbVertices, params);
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


/* 
 * Convert polynomial from long long C-style matrix to Ginac expression
 *
 *	m: polynomial matrix
 *	nbRows, nbColumns: number of rows, columns
 */
int polyConvert(long long *m, unsigned int nbRows, unsigned int nbColumns) 
{
	unsigned int nbVariables = nbColumns - 2;
	ex p;

	// setting global variables
	Vars = constructVariableVector(nbVariables, "v");
       
	p = convertPolynomial(m, nbRows, nbColumns, Vars);
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
int polyConvertParameters(long long *m, unsigned int nbRows, unsigned int nbColumns,
			  long long **llPolynomialCoefficients, unsigned int *llRows, unsigned int *llColumns,
			  const exvector& params)
{
	unsigned int nbVariables = nbColumns;
	unsigned int nbCoefficients = nbRows;
	ex p;

	// setting global variables
	Vars = constructVariableVector(nbVariables, "v");
       
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
				t += pow(Vars[j], (long int) m[i*nbColumns+j]);
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


/* Reads the polynomial matrix, converts it to long long precision and calls ginac functions */
Matrix *readPolynomial(unsigned int nbVariables, const exvector& params)
{
	Matrix *polynomial;

	polynomial = Matrix_Read();

#ifdef DEBUG
	/* Print the polynomial matrix */
	printf("================================\n");
	printf("Polynomial: \n");
	Matrix_Print(stdout, P_VALUE_FMT, polynomial);
	printf("================================\n");
#endif

	long long *matrix = matrix2longlong(polynomial);
	// parameters in the polynomial coefficients
	if(polynomial->NbColumns == nbVariables) {
		unsigned int nbCoefficients = polynomial->NbRows;;
		unsigned int i;

		// FIXME: free
		long long **llPolynomialCoefficients = (long long **) calloc(sizeof(long *), nbCoefficients);
		Matrix **mPolynomialCoefficients = (Matrix **) calloc(sizeof(Matrix *), nbCoefficients);
		unsigned int *llRows = (unsigned int *) calloc(sizeof(unsigned int), nbCoefficients);
		unsigned int *llColumns = (unsigned int *) calloc(sizeof(unsigned int), nbCoefficients);;

		for(i = 0; i < nbCoefficients; i++) {
			// read the matrix and set rows and columns number.
			mPolynomialCoefficients[i] = Matrix_Read();
			llRows[i] = mPolynomialCoefficients[i]->NbRows;
			llColumns[i] = mPolynomialCoefficients[i]->NbColumns;
#ifdef DEBUG
			/* Print the i coefficient matrix */
			printf("================================\n");
			printf("Coefficient i: \n");
			Matrix_Print(stdout, P_VALUE_FMT, mPolynomialCoefficients[i]);
			printf("================================\n");
#endif
			llPolynomialCoefficients[i] = matrix2longlong(mPolynomialCoefficients[i]);
			Matrix_Free(mPolynomialCoefficients[i]);
		}
		free(mPolynomialCoefficients);
		polyConvertParameters(matrix, polynomial->NbRows, polynomial->NbColumns, llPolynomialCoefficients
				      , llRows, llColumns, params);
		for(i = 0; i < nbCoefficients; i++) {
			free(llPolynomialCoefficients[i]);
		}
		free(llRows);
		free(llColumns);
		free(llPolynomialCoefficients);

	} else {
		polyConvert(matrix, polynomial->NbRows, polynomial->NbColumns);
	}
	free(matrix);

	return polynomial;
}
