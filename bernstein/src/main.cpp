#include <stdio.h>
#include <stdlib.h>

#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"

#include <bernstein/bernstein.h>
#include <bernstein/maximize.h>
#include "polynomial.h"

#define MAXRAYS 1000

using namespace std;
using namespace GiNaC;
using namespace bernstein;

static ex readPolynomial(const exvector& vars, const exvector& params);
static void printCoefficients(lst coeffs);
static int printMaxMinCoefficient(Polyhedron *VD, lst coeffs, const exvector& Params);

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
	param_name = Read_ParamNames(stdin, nb_var+nb_param);
	vars = constructParameterVector(param_name, nb_var);
	params = constructParameterVector(param_name+nb_var, nb_param);

	polynomial = readPolynomial(vars, params);

	Matrix_Free(a);
	Matrix_Free(b);

	/* Find the parametrized domains */
	PP = Polyhedron2Param_Domain(A,B,MAXRAYS);
	for(Q=PP->D;Q;Q=Q->next) {
		Polyhedron *VD;
		lst coeffs;

		printf("\nDomain: \n");
		VD = DomainSimplify(Q->Domain, B, MAXRAYS);
		Print_Domain(stdout, VD, param_name+nb_var);
		matrix VM = domainVertices(PP, Q, params);
		coeffs = bernsteinExpansion(VM, polynomial, vars, params);
		printCoefficients(coeffs);
		printMaxMinCoefficient(Q->Domain, coeffs, params);
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
	char buffer[1024], *s;
	lst allvars;
	ex p;

	for (int i = 0; i < vars.size(); ++i)
		allvars.append(vars[i]);
	for (int i = 0; i < params.size(); ++i)
		allvars.append(params[i]);

	do 
		s = fgets(buffer, 1024, stdin);
	while (s && (s[0] == '#' || s[0] == '\n'));

	if (!s)
		return 0;

	try {
		p = ex(string(s), allvars);
	} catch (exception &p) {
		cerr << p.what() << endl;
		return 0;
	}

	return p;
}


int printMaxMinCoefficient(Polyhedron *VD, lst coeffs, const exvector& Params)
{
	cout << "\tMinimum coefficient(s): " << minimize(VD, coeffs, Params) << endl;
	cout << "\tMaximum coefficient(s): " << maximize(VD, coeffs, Params) << endl;
	return 0;
}
