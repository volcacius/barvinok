#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include <polylib/polylibgmp.h>
#undef value_compare
#undef divide

#include "bernstein.h"
#include "bernstein++.h"

using namespace GiNaC;

// global to be used on the coefficient max
Polyhedron *VD = NULL;	// current validity domain

/* main function */
int main(void) {
	Matrix *a, *b, *polynomial;
	Polyhedron *A, *B;		// initial matrices
	char **param_name;	// name of the parameters
	exvector params;

	Param_Polyhedron *PP;
	Param_Domain   *Q;

	unsigned int nb_param;

	printf("\n===============================================\n");

	a = Matrix_Read();
	A = Constraints2Polyhedron(a,200);

	b = Matrix_Read();
	B = Constraints2Polyhedron(b,200);

	/* Set the number of parameters */
	nb_param	= B->Dimension;

	/* Read the name of the parameters */
	param_name = Read_ParamNames(stdin, nb_param);
	params = constructParameterVector(param_name, nb_param);

	polynomial = readPolynomial(a->NbColumns-b->NbColumns, params);

	Matrix_Free(a);
	Matrix_Free(b);

	/* Find the parametrized domains */
	PP = Polyhedron2Param_SimplifiedDomain(&A,B,MAXRAYS,NULL,NULL);
	for(Q=PP->D;Q;Q=Q->next) {
		printDomain(Q, B, param_name);
		doExpansion(PP, Q, params);
		printf("\n\n===============================================\n");
	}

	Domain_Free(A);
	Domain_Free(B);
	Param_Polyhedron_Free(PP);
	free(param_name);

	Matrix_Free(polynomial);

	Domain_Free(VD);

	return 0;
} /* main */
