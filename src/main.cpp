#include <stdio.h>
#include <stdlib.h>

#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"

#include "bernstein.h"
#include "bernstein-expansion.h"

using namespace std;
using namespace GiNaC;

static ex readPolynomial(const exvector& vars, const exvector& params);
static void printCoefficients(lst coeffs);
static long long value2longlong(Value p);
static long long *matrix2longlong(Matrix *M);
static ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns,
			    const exvector& params);
static ex polyConvertParameters(long long *m, unsigned int nbRows, 
			  unsigned int nbColumns, 
			  long long **llPolynomialCoefficients, unsigned int *llRows, 
			  unsigned int *llColumns, const exvector& vars,
			  const exvector& params);

/* main function */
int main(void) {
	Matrix *a, *b;
	Polyhedron *A, *B;		// initial matrices
	char **param_name;	// name of the parameters
	exvector params, vars;
	ex polynomial;

	Param_Polyhedron *PP;
	Param_Domain   *Q;

	unsigned int nb_param, nb_var;

	printf("\n===============================================\n");

	a = Matrix_Read();
	A = Constraints2Polyhedron(a,200);

	b = Matrix_Read();
	B = Constraints2Polyhedron(b,200);

	/* Set the number of parameters */
	nb_param	= B->Dimension;
	nb_var		= A->Dimension - B->Dimension;

	/* Read the name of the parameters */
	param_name = Read_ParamNames(stdin, nb_param);
	params = constructParameterVector(param_name, nb_param);
	vars = constructVariableVector(nb_var, "v");

	polynomial = readPolynomial(vars, params);

	Matrix_Free(a);
	Matrix_Free(b);

	/* Find the parametrized domains */
	PP = Polyhedron2Param_SimplifiedDomain(&A,B,MAXRAYS,NULL,NULL);
	for(Q=PP->D;Q;Q=Q->next) {
		Polyhedron *VD;
		lst coeffs;

		printf("\nDomain: \n");
		VD = DomainSimplify(Q->Domain, B, MAXRAYS);
		Print_Domain(stdout, VD, param_name);
		matrix VM = domainVertices(PP, Q, params);
		coeffs = bernsteinExpansion(VM, polynomial, vars, params);
		printCoefficients(coeffs);
		getMaxMinCoefficient(Q->Domain, coeffs, params);
		Domain_Free(VD);
		printf("\n\n===============================================\n");
	}

	Domain_Free(A);
	Domain_Free(B);
	Param_Polyhedron_Free(PP);
	free(param_name);

	return 0;
} /* main */


void printCoefficients(lst coeffs)
{
	cout << "-----------------------------------------------" << endl;
	cout << "Coefficients: " << endl << endl;

	for (lst::const_iterator i = coeffs.begin(); i != coeffs.end(); ++i)
	    cout << "\t" << *i << endl;
}


/* Reads the polynomial matrix, converts it to long long precision and calls ginac functions */
ex readPolynomial(const exvector& vars, const exvector& params)
{
	ex p;
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
	if(polynomial->NbColumns == vars.size()) {
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
		p = polyConvertParameters(matrix, polynomial->NbRows, 
					  polynomial->NbColumns,
					  llPolynomialCoefficients,
					  llRows, llColumns, vars, params);
		for(i = 0; i < nbCoefficients; i++) {
			free(llPolynomialCoefficients[i]);
		}
		free(llRows);
		free(llColumns);
		free(llPolynomialCoefficients);

	} else {
		p = convertPolynomial(matrix, polynomial->NbRows, 
				      polynomial->NbColumns, vars);
	}
	delete [] matrix;
	Matrix_Free(polynomial);

	return p;
}


/* Converts a value (p) to a long long integer */
long long value2longlong(Value p)
{
	char *str; 
	long long l;

	str = mpz_get_str(0,10, p);
	l = atoll(str);
	free(str);

	return l;
}


/* Converts a Matrix from the polylib format to long long* format */
long long *matrix2longlong(Matrix *M)
{
 	long long *matrix = new (long long)[M->NbRows * M->NbColumns];
	unsigned int nr,nc;
	unsigned int j,k;
	Value *p;

	p=*(M->p);
	nr = M->NbRows;
	nc = M->NbColumns;

	/* Copy matrix to long long */
	for(j = 0; j < nc; j++) {
		for(k = 0; k < nr; k++) {
			// TODO: loosing precision to long long
			long long l = value2longlong(*p++);
			matrix[nr*j+k] = l;
		}
	}
	return matrix;
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
