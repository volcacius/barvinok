#include <stdio.h>
#include <stdlib.h>

#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"

#include <bernstein/bernstein.h>
#include "polynomial.h"

#define MAXRAYS 1000

using namespace std;
using namespace GiNaC;
using namespace bernstein;

static ex readPolynomial(const exvector& vars, const exvector& params);
static void printCoefficients(lst coeffs);
static unsigned int findMaxDegree(lst polynomials, ex var);
static bool linearCoefficients(lst coeffs, const exvector &Params);
static bool generateMaxConstraints(Polyhedron *VD, lst coeffs, 
				   const exvector &Params, unsigned int max);
static bool generateMinConstraints(Polyhedron *VD, lst coeffs, 
				   const exvector &Params, unsigned int min);
static bool generatePositiveNegativeConstraints(Polyhedron *VD, bool positive);
static void ex2longlongRow(long long *M, unsigned int row, polynomial &difference,
			   const exvector& Params);
static ex constantTerm(ex poly, const exvector &Params);
static ex calculateLCM(ex poly, const exvector &Params, ex &cte);
static ex getMaxMinCoefficient1Param(lst coeffs, unsigned int maxDegree, 
				     ex Param, bool max, bool positive);
static bool calculateDirection(bool max, bool positive, bool even);
static int getMaxMinCoefficient(Polyhedron *VD, lst coeffs, const exvector& Params);
static unsigned  checkConstraint(Polyhedron *VD, long long *M, 
			  unsigned int rows, unsigned int columns);
static Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns);

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


int getMaxMinCoefficient(Polyhedron *VD, lst coeffs, const exvector& Params)
{
	if(Params.size() == 1) {
		unsigned maxDegree = findMaxDegree(coeffs, Params[0]);
		// check if the parameter is positive
		if(generatePositiveNegativeConstraints(VD, true)) {
			cout << endl << Params[0] << " >= 0" << endl;
			ex m = getMaxMinCoefficient1Param(coeffs, maxDegree, Params[0], true, true);
			cout << "\tMaximum coefficient: " << m << endl;
			m = getMaxMinCoefficient1Param(coeffs, maxDegree, Params[0], false, true);
			cout << "\tMinimum coefficient: " << m << endl;
		}
		// check if the parameter is negative
		if(generatePositiveNegativeConstraints(VD, false)) {
			cout << endl << Params[0] << " < 0" << endl;
			ex m = getMaxMinCoefficient1Param(coeffs, maxDegree, Params[0], true, false);
			cout << "\tMaximum coefficient: " << m << endl;
			m = getMaxMinCoefficient1Param(coeffs, maxDegree, Params[0], false, false);
			cout << "\tMinimum coefficient: " << m << endl;
		}
	}
	if(Params.size() > 1 && linearCoefficients(coeffs, Params)) {
#ifdef DEBUG
		cout << "=================================================" << endl;
		cout << "Linear coefficients:" << endl << endl;
#endif
		for(unsigned int k = 0; k < coeffs.nops(); k++) {
#ifdef DEBUG
			cout << "#################################################" << endl;
			cout << "Proposing max: " << coeffs[k] << endl << endl;
#endif
			if(generateMaxConstraints(VD, coeffs, Params, k)) {

				cout << "\tMaximum coefficient: " << coeffs[k] << endl;
			}
#ifdef DEBUG
			cout << "#################################################" << endl;
			cout << "Proposing min: " << coeffs[k] << endl << endl;
#endif

			if(generateMinConstraints(VD, coeffs, Params, k)) {
				cout << "\tMinimum coefficient: " << coeffs[k] << endl;
			}

		}
	}

	return 0;
}

/*
 * Check if all the coefficients are linear on the parameters
 *
 *	coeffs: coefficients list
 *	Params: parameters matrix
 *	nbParams: number of parameters
 */
bool linearCoefficients(lst coeffs, const exvector &Params)
{
	for (size_t i = 0; i < coeffs.nops(); ++i)
		for (unsigned int j = 0; j < Params.size(); j++)
			if (coeffs[i].degree(Params[j]) > 1)
				return false;
	return true;
}



/*
 * Generate constraints where max coeff is greater than the others and 
 * call function that checks if it holds
 *
 *	coeffs: coefficients list
 *	Params: parameters matrix
 *	nbParams: number of parameters
 *	max: proposed maximum
 */
bool generateMaxConstraints(Polyhedron *VD, lst coeffs, const exvector &Params, 
			    unsigned int max)
{
	lst differences;

	long long *M; // diferences matrix
	M = (long long *) calloc(sizeof(long long), 
				 (coeffs.nops()-1) * (Params.size()+2));		

	unsigned int row = 0;
	for (size_t i = 0; i < coeffs.nops(); ++i) {
		if(i != max) {
			polynomial difference(coeffs[max] - coeffs[i]);
#ifdef DEBUG
			cout << "Diff(0): " << difference.term(0) << endl;
			cout << "C_max - C_i -1: " << difference << endl << endl;
#endif
			ex2longlongRow(M, row, difference, Params);
			row++;
		}
	}
	bool retval = checkConstraint(VD, M, coeffs.nops()-1, Params.size()+2);
	free(M);
	return retval;
}


/*
 * Generate constraints where min coeff is greater than the others and 
 * call function that checks if it holds
 *
 *	coeffs: coefficients list
 *	Params: parameters matrix
 *	nbParams: number of parameters
 *	max: proposed maximum
 */
bool generateMinConstraints(Polyhedron *VD, lst coeffs, const exvector &Params, 
			    unsigned int min)
{
	lst differences;

	long long *M; // diferences matrix
	M = (long long *) calloc(sizeof(long long), 
				 (coeffs.nops()-1) * (Params.size()+2));		

	unsigned int row = 0;
	for (size_t i = 0; i < coeffs.nops(); ++i) {
		if(i != min) {
			polynomial difference(coeffs[i] - coeffs[min]);
#ifdef DEBUG
			cout << "C_max - C_i -1: " << difference << endl << endl;
#endif
			ex2longlongRow(M, row, difference, Params);
			row++;
		}
	}
	bool retval = checkConstraint(VD, M, coeffs.nops()-1, Params.size()+2);
	free(M);
	return retval;
}


/*
 * Generate (and check) constraints where the parameter is positive/negative
 *
 *	coeffs: coefficients list
 */
bool generatePositiveNegativeConstraints(Polyhedron *VD, bool positive)
{
	long long *M;
	M = (long long *) calloc(sizeof(long long), 1 * 3);

	M[0] = 1;
	if(positive) {
		M[1] = 1;
		M[2] = 0;
	} else {
		M[1] = -1;
		M[2] = -1;
	}

	bool retval = checkConstraint(VD, M, 1, 3);
	free(M);
	return retval;
}



/*
 * Convert ex to polylib-like matrix format row
 *
 *	M: polylib-like matrix format
 *	row: row to fill
 *	difference: polynome
 *	Params: parameters matrix
 *	nbParams: number of parameters
 */
void ex2longlongRow(long long *M, unsigned int row, polynomial &difference, 
		    const exvector& Params)
{
	unsigned int columns = Params.size()+2;
	ex cte = constantTerm(difference, Params);
	ex lcm = calculateLCM(difference, Params, cte);

	M[row*columns] = 1;		// greater or equal

	// TODO: losing precision
	long longCte = (ex_to<numeric>(cte).to_long());
	M[row*columns+columns-1] = longCte;	// cte

	// difference with integer coefficients
	polynomial intDifference(difference * lcm);

	// now, fill the matrix row
	for (size_t i = 0; i < intDifference.nbTerms(); ++i) {
		for(unsigned int j = 0; j < Params.size(); j++) {
			if(intDifference.term(i).degree(Params[j]) == 1) {
#ifdef DEBUG
				cout << "difference.term(" << i << "): " << (intDifference.term(i)).coeff(Params[j]) << endl;
#endif
				// TODO: losing precision
				long val = (ex_to<numeric>(intDifference.term(i).coeff(Params[j]))).to_long();
#ifdef DEBUG
				cout << "(i,j): " << row << "," << j+1 << ": " << val << endl;
#endif
				M[row*columns+j+1] = val;
			}
		}
	}
}


/*
 * Find the constant term of a multi-parameter polynome
 *
 *	p: polynome
 *	Params: parameters matrix
 *	nbParams: number of parameters
 */
ex constantTerm(ex p, const exvector &Params)
{
	ex constants;
	polynomial poly(p);

#ifdef DEBUG
	cout << "Input poly: " << poly << endl;
#endif

	for (size_t i = 0; i < poly.nbTerms(); ++i) {
#ifdef DEBUG
		cout << "poly(" << i << "): " << poly.term(i) << endl;
#endif
		bool cnt = true;	// is this term constant?

		for(unsigned int j = 0; j < Params.size(); j++) {
			if(poly.term(i).degree(Params[j]) == 1) {
				cnt = false;
			}
		}
		if(cnt) {
			constants += poly.term(i);
		}
	}
#ifdef DEBUG
	cout << "Constants: " << constants << endl;
#endif
	return constants;
}



/*
 * Find the lcm between the coefficients denominators
 *
 *	p: polynome
 *	Params: parameters matrix
 *	nbParams: number of parameters
 *	cte: constant term
 */
ex calculateLCM(ex p, const exvector &Params, ex &cte)
{
	ex retval = 1;
	polynomial poly(p);

	for (size_t i = 0; i < poly.nbTerms(); ++i) {
		for(unsigned int j = 0; j < Params.size(); j++) {
			if(poly.term(i).degree(Params[j]) == 1) {
				retval = lcm(retval, poly.coeff(Params[j]).denom());
			}
		}
	}
	if(cte != 0) {
		retval = lcm(retval, cte.denom());
	}
#ifdef DEBUG
	cout << "LCM: " << retval << endl;
#endif
	return retval;
}

/*
 * Get the maximum or minimum coefficient in the one parameter case
 *
 *	coeffs: coefficients list
 *	maxDegree: maximum degree of the polynomial
 *	Params: parameters matrix
 *	max: flag that indicates if we are searching for max or min
 *	positive: flag that indicates if the parameter is negative or positive
 */
ex getMaxMinCoefficient1Param(lst coeffs, unsigned int maxDegree, ex Param, bool max, bool positive)
{
	bool reverse;		// indicates the direction of the coeff comparation
	int d = maxDegree;
	bool even = (d%2==0);
	lst m;
	m.append(coeffs[0]);

#ifdef DEBUG
	cout << "Max/Min: " << max << endl;
	cout << "Positive/Negative: " << positive << endl;
	cout << "Even/Odd: " << even << endl;

	cout << "Deg: " << d << endl;
	cout << "Input Coeffs: " << coeffs <<endl;
	cout << "Maxs/Mins: " << m <<endl;

#endif
	for (size_t i = 1; i < coeffs.nops(); ++i) {
		ex degreeDCoeff = coeffs[i].coeff(Param, d);
#ifdef DEBUG
		cout << "Coef: " << coeffs[i] << " -> " << degreeDCoeff << endl;
#endif
		reverse = calculateDirection(max, positive, even);
		if(!reverse) {
			if(degreeDCoeff >= m[0].coeff(Param, d)) {
				if(degreeDCoeff > m[0].coeff(Param, d)) {
					m.remove_all();
				}
				m.append(coeffs[i]);
			}
		} else {
			if(degreeDCoeff <= m[0].coeff(Param, d)) {
				if(degreeDCoeff < m[0].coeff(Param, d)) {
					m.remove_all();
				}
				m.append(coeffs[i]);
			}
		}

#ifdef DEBUG
		cout << "Maxs/Mins: " << m <<endl;
#endif
	}
	if(m.nops() > 1 && d > 0) {
		return getMaxMinCoefficient1Param(m, d-1, Param, max, positive);
	} else {
		return m[0];
	}
}

/*
 * Calculate the direction of the comparation in the one parameter case
 *
 *	max: flag that indicates if we are searching for max or min
 *	positive: flag that indicates if the parameter is negative or positive
 *	even: flag that indicates if the actual degree is even or odd
 */
bool calculateDirection(bool max, bool positive, bool even)
{
	bool reverse;

	if(positive) {
		if(max) {
			reverse = false;
		} else {
			reverse = true;
		}
	}
	// negative
	if(!positive) {
		// odd
		if(!even) {
			if(max) {
				reverse = true;
			} else {
				reverse = false;
			}
		}
		if(even) {
			if(max) {
					reverse = false;
			} else {
				reverse = true;
			}
		}
	}
	return reverse;
}


unsigned int findMaxDegree(lst polylst, ex var)
{
	unsigned max = 0;
	for (lst::const_iterator i = polylst.begin(); i != polylst.end(); ++i) {
		unsigned degree = i->degree(var);
		if (degree > max)
		    max = degree;
	}
	return max;
}


/* Converts a *longlong matrix to polylib format */
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns)
{
	Matrix *retval;
	unsigned int i, j;

	retval = Matrix_Alloc(rows, columns);
	for(i = 0; i < rows; i++) {
		for(j = 0; j < columns; j++) {
			value_set_si(retval->p[i][j], (int) M[i*columns+j]);
		}
	}
	return retval;
}


/* Check  if a given set of constraints (M matrix) holds or not in the validity domain */
unsigned  checkConstraint(Polyhedron *VD, long long *M, 
			  unsigned int rows, unsigned int columns)
{
	Matrix *maxConstraints;
	Polyhedron *mC, *newB;

	maxConstraints = longlong2polylib(M, rows, columns);

#ifdef DEBUG
	printf("Max/Min Constraints (polylib Format): \n");
	Matrix_Print(stdout, P_VALUE_FMT, maxConstraints);
	printf("\n");

	printf("Original Validity Domain: \n");
	Print_Domain(stdout, VD, param_name);
	printf("\n");
#endif

	mC =Constraints2Polyhedron(maxConstraints, MAXRAYS);

	if(!PolyhedronIncludes(mC, VD)) {
#ifdef DEBUG
		printf("Proposed Coefficient is not the max/min/pos/neg.\n");
#endif
		return 0;
	} else {
#ifdef DEBUG
		printf("Proposed Coefficient is the max/min/pos/neg.\n");
#endif
		return 1;
	}

	/* free */
	Matrix_Free(maxConstraints);
	Domain_Free(mC);
	Domain_Free(newB);

	return 0;
}
