#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include <polylib/polylibgmp.h>
#undef value_compare
#undef divide

#include "bernstein.h"
#include "bernstein++.h"

using namespace GiNaC;

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

	polynomial = readPolynomial(nb_var, vars, params);

	Matrix_Free(a);
	Matrix_Free(b);

	/* Find the parametrized domains */
	PP = Polyhedron2Param_SimplifiedDomain(&A,B,MAXRAYS,NULL,NULL);
	for(Q=PP->D;Q;Q=Q->next) {
		Polyhedron *VD;
		printf("\nDomain: \n");
		VD = DomainSimplify(Q->Domain, B, MAXRAYS);
		Print_Domain(stdout, VD, param_name);
		doExpansion(PP, Q, polynomial, vars, params);
		Domain_Free(VD);
		printf("\n\n===============================================\n");
	}

	Domain_Free(A);
	Domain_Free(B);
	Param_Polyhedron_Free(PP);
	free(param_name);

	return 0;
} /* main */
