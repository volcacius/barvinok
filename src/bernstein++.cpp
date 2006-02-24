/* 
 *	Bernstein Expansion
 *
 *	- c to c++ functions
 */



#include <iostream>
#include <string>
#include <cln/cln.h>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

#include "bernstein++.h"
#include "bernstein.h"
#include "bernstein-expansion.h"


static numeric value2numeric(Value v)
{
    int sa = v[0]._mp_size;
    if (!sa)
	return 0;
    int abs_sa = sa < 0 ? -sa : sa;
    cln::cl_I res = 0;
    for (int i = abs_sa-1; i >= 0; --i) {
	res = res << GMP_LIMB_BITS;
	res = res + v[0]._mp_d[i];
    }
    return numeric(sa < 0 ? -res : res);
}


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
lst doExpansion(Param_Polyhedron *PP, Param_Domain *Q, ex polynomial, 
		const exvector& vars, const exvector& params)
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
			for (unsigned j = 0; j < nbColumns-2; j++)
				t += value2numeric(V->Vertex->p[i][j]) * params[j];
			t += value2numeric(V->Vertex->p[i][nbColumns-2]);
			t /= value2numeric(V->Vertex->p[i][nbColumns-1]);
#ifdef DEBUG
			cout << "T: " << t << endl;
#endif
			VM(v, i) = t;
		}
		++v;
	END_FORALL_PVertex_in_ParamPolyhedron;

	return bernsteinExpansion(VM, polynomial, vars, params);
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
