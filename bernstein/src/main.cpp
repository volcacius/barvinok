#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"

#include <bernstein/bernstein.h>
#include <bernstein/maximize.h>
#include <bernstein/piecewise_lst.h>

#define MAXRAYS 1000

using namespace std;
using namespace GiNaC;
using namespace bernstein;

static ex readPolynomial(const exvector& vars, const exvector& params);
static void printCoefficients(lst coeffs);
static int printMaxMinCoefficient(Polyhedron *VD, lst coeffs, const exvector& Params);

static char *readLine(void)
{
	static char buffer[1024];
	char *s;

	do 
		s = fgets(buffer, 1024, stdin);
	while (s && (s[0] == '#' || s[0] == '\n'));

	return s;
}

static void readExpected(const exvector& vars, const exvector& params,
			 piecewise_lst **exp_all,
			 piecewise_lst **exp_min,
			 piecewise_lst **exp_max)
{
	int n = 0;
	char *s;
	piecewise_lst *pl_all, *pl_min, *pl_max;

	pl_all = new piecewise_lst(params);
	pl_min = new piecewise_lst(params);
	pl_max = new piecewise_lst(params);

	s = readLine();
	assert(s);

	sscanf(s, "%d", &n);

	for (int i = 0 ; i < n ; ++i) {
		ex all, min, max;
		Polyhedron *D;
		Matrix *M = Matrix_Read();

		D = Constraints2Polyhedron(M, MAXRAYS);
		Matrix_Free(M);

		all = readPolynomial(vars, params);
		min = readPolynomial(vars, params);
		max = readPolynomial(vars, params);
		assert(is_a<lst>(all));
		assert(is_a<lst>(min));
		assert(is_a<lst>(max));
		pl_all->add_guarded_lst(Polyhedron_Copy(D), ex_to<lst>(all));
		pl_min->add_guarded_lst(Polyhedron_Copy(D), ex_to<lst>(min));
		pl_max->add_guarded_lst(D, ex_to<lst>(max));
	}
	*exp_all = pl_all;
	*exp_min = pl_min;
	*exp_max = pl_max;
}

/* main function */
int main(int argc, char *argv[])
{
	Matrix *a, *b;
	Polyhedron *A, *B;		// initial matrices
	char **param_name;	// name of the parameters
	exvector params, vars;
	ex polynomial;
	piecewise_lst *exp_all, *exp_min, *exp_max;
	piecewise_lst *got_all, *got_min, *got_max;

	Param_Polyhedron *PP;
	Param_Domain   *Q;

	unsigned int nb_param, nb_var;

	int c;
	int verify = 0;

	while ((c = getopt(argc, argv, "T")) != -1) {
		switch(c) {
		case 'T':
			verify = 1;
			break;
		}
	}

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

	if (verify) {
		readExpected(vars, params, &exp_all, &exp_min, &exp_max);
		got_all = new piecewise_lst(params);
		got_min = new piecewise_lst(params);
		got_max = new piecewise_lst(params);
	}

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

		if (!verify)
			continue;

		got_all->add_guarded_lst(Polyhedron_Copy(Q->Domain), coeffs);
		got_min->add_guarded_lst(Polyhedron_Copy(Q->Domain),
					 minimize(Q->Domain, coeffs, params));
		got_max->add_guarded_lst(Polyhedron_Copy(Q->Domain),
					 maximize(Q->Domain, coeffs, params));
	}

	Domain_Free(A);
	Domain_Free(B);
	Param_Polyhedron_Free(PP);
	free(param_name);

	if (verify) {
		if (!got_all->is_equal(*exp_all)) {
			cerr << "expected:" << endl;
			cerr << *exp_all << endl;
			cerr << "got:" << endl;
			cerr << *got_all << endl;
			return 1;
		}
		if (!got_min->is_equal(*exp_min)) {
			cerr << "expected:" << endl;
			cerr << *exp_min << endl;
			cerr << "got:" << endl;
			cerr << *got_min << endl;
			return 1;
		}
		if (!got_max->is_equal(*exp_max)) {
			cerr << "expected:" << endl;
			cerr << *exp_max << endl;
			cerr << "got:" << endl;
			cerr << *got_max << endl;
			return 1;
		}

		delete exp_all;
		delete exp_max;
		delete exp_min;
		delete got_all;
		delete got_max;
		delete got_min;
	}

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
	char *s;
	lst allvars;
	ex p;

	for (int i = 0; i < vars.size(); ++i)
		allvars.append(vars[i]);
	for (int i = 0; i < params.size(); ++i)
		allvars.append(params[i]);

	s = readLine();
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
