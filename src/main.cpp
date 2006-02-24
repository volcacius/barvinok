#include <stdio.h>
#include <stdlib.h>

#include <ginac/ginac.h>
#include <gmp.h>

#define polynomial polylib_polynomial
#include <polylib/polylibgmp.h>
#undef polynomial
#undef value_compare
#undef divide

using namespace GiNaC;

#include "bernstein.h"
#include "bernstein++.h"
#include "bernstein-expansion.h"

static ex readPolynomial(const exvector& vars, const exvector& params);

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
		coeffs = doExpansion(PP, Q, polynomial, vars, params);
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
	free(matrix);
	Matrix_Free(polynomial);

	return p;
}
