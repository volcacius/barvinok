/* 
 *	Bernstein Expansion
 *
 *	- ginac functions
 */


#include <iostream>

#include <string>

#include "bernstein-expansion.h"

static unsigned int findMaxDegree(ex polynomial, const exvector& Vars);

/*
 * Do the Bernstein Expansion 
 *
 *	P: parameters matrix
 *	poly: polynomial expression
 *	V: variables matrix
 *	nbVar: number of variables
 *	nbVert: number of vertices
 *	maxDegree: max multi-degree of the polynomial
 */
int bernsteinExpansion(matrix &P, ex &poly, const exvector& V,
		       unsigned int nbVert, const exvector& Params)
{
	unsigned maxDegree = findMaxDegree(poly, V);
	matrix A = getAiMatrix(nbVert);

#ifdef DEBUG
	cout << endl << "Polynomial: " << poly << endl << endl;
	cout << P << endl << endl;
	cout << A << endl << endl;
#endif

	// obtain the variables value on the basis and replace
	ex variables = evalm(A * P);
	ex polynom = replaceVariablesInPolynomial(poly, V, variables);

#ifdef DEBUG
	cout << variables << endl << endl;
	cout << "Preliminar Expansion: " << polynom << endl << endl;
#endif

	ex basis = getBasis(nbVert, A);
	ex basisPowered = pow(basis, maxDegree);
	polynomial expandedBasis = basisPowered.expand();

#ifdef DEBUG
	cout << "Basis: " << basis << endl<< endl;
	cout << "Powered Basis: " << basisPowered << endl;
	cout << "Expanded Basis: " << expandedBasis << endl << endl;
#endif

	// monomials to n degree
	polynomial p(polynom);
	ex maxDegreePolynomial = powerMonomials(p, A, nbVert
						, maxDegree, basis);

	// get the coefficients
	lst coeffs = getCoefficients(maxDegreePolynomial, expandedBasis, nbVert, A);
	if(Params.size() == 1) {
		// check if the parameter is positive
		if(generatePositiveNegativeConstraints(true)) {
			cout << endl << Params[0] << " >= 0" << endl;
			ex m = getMaxMinCoefficient1Param(coeffs, maxDegree, Params[0], true, true);
			cout << "\tMaximum coefficient: " << m << endl;
			m = getMaxMinCoefficient1Param(coeffs, maxDegree, Params[0], false, true);
			cout << "\tMinimum coefficient: " << m << endl;
		}
		// check if the parameter is negative
		if(generatePositiveNegativeConstraints(false)) {
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
			if(generateMaxConstraints(coeffs, Params, k)) {

				cout << "\tMaximum coefficient: " << coeffs[k] << endl;
			}
#ifdef DEBUG
			cout << "#################################################" << endl;
			cout << "Proposing min: " << coeffs[k] << endl << endl;
#endif

			if(generateMinConstraints(coeffs, Params, k)) {
				cout << "\tMinimum coefficient: " << coeffs[k] << endl;
			}

		}
	}

	return 0;
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
				cout << "difference.term(" << i << "): " << (intDifference.term(i)).coeff(Params(0,j)) << endl;
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
 * Generate constraints where max coeff is greater than the others and 
 * call function that checks if it holds
 *
 *	coeffs: coefficients list
 *	Params: parameters matrix
 *	nbParams: number of parameters
 *	max: proposed maximum
 */
bool generateMaxConstraints(lst coeffs, const exvector &Params, unsigned int max)
{
	lst differences;

	long long *M; // diferences matrix
	M = (long long *) calloc(sizeof(long long), 
				 (coeffs.nops()-1) * (Params.size()+2));		

	unsigned int row = 0;
	for (size_t i = 0; i < coeffs.nops(); ++i) {
		if(i != max) {
			polynomial difference(coeffs[max] - coeffs[i] -1);
#ifdef DEBUG
			cout << "Diff(0): " << difference.term(0) << endl;
			cout << "C_max - C_i -1: " << difference << endl << endl;
#endif
			ex2longlongRow(M, row, difference, Params);
			row++;
		}
	}
	bool retval = checkConstraint(M, coeffs.nops()-1, Params.size()+2);
	free(M);
	return retval;
}

/*
 * Generate (and check) constraints where the parameter is positive/negative
 *
 *	coeffs: coefficients list
 */
bool generatePositiveNegativeConstraints(bool positive)
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

	bool retval = checkConstraint(M, 1, 3);
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
bool generateMinConstraints(lst coeffs, const exvector &Params, unsigned int min)
{
	lst differences;

	long long *M; // diferences matrix
	M = (long long *) calloc(sizeof(long long), 
				 (coeffs.nops()-1) * (Params.size()+2));		

	unsigned int row = 0;
	for (size_t i = 0; i < coeffs.nops(); ++i) {
		if(i != min) {
			polynomial difference(coeffs[i] - coeffs[min] -1);
#ifdef DEBUG
			cout << "C_max - C_i -1: " << difference << endl << endl;
#endif
			ex2longlongRow(M, row, difference, Params);
			row++;
		}
	}
	bool retval = checkConstraint(M, coeffs.nops()-1, Params.size()+2);
	free(M);
	return retval;
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
	bool retval = true;

	for (size_t i = 0; i < coeffs.nops(); ++i) {
		for(unsigned int j = 0; j < Params.size(); j++) {
			retval = retval && (coeffs[i].degree(Params[j]) <= 1);
		}
	}
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


/*
 * Construct A_i matrix
 *
 *	nbVert: number of vertices
 */
matrix getAiMatrix(unsigned int nbVert)
{
	matrix A(1, nbVert); 		// a_i matrix
	symbol *pSym;		// pointer to a_i
	for(unsigned int i = 0; i < nbVert; i++) {
		pSym = new symbol("a" + int2String(i));
#ifdef DEBUG
		cout << "A: " << *pSym << endl;
#endif
		A(0,i) = *pSym;
		delete pSym;
	}
	return A;
}


/*
 * Construct the basis
 *
 *	A: a_i matrix
 *	nbVert: number of vertices
 */
ex getBasis(unsigned int nbVert, matrix &A)
{
	ex basis;
	for(unsigned int i = 0; i < nbVert; i++) {
		basis += A(0,i);
	}
	return basis;
}


/*
 * Finds the less than maxDegree monomials and multiply them
 * by the basis
 *
 *	polynomial: original polynomial
 *	A: a_i matrix
 *	nbVert: number of vertices
 *	maxDegree: maximal polynomial multidegree
 *	basis: basis of the polytope
 */
ex powerMonomials(polynomial &poly, matrix &A, unsigned int nbVert
		  , unsigned int maxDegree, ex &basis)
{
	ex maxDegreePolynomial;
#ifdef DEBUG
	cout << "- Degree --------------------------------------" << endl;
#endif
	for (size_t i = 0; i != poly.nbTerms(); ++i) {
		unsigned int degree = 0;		// degree of the monomial

		for(unsigned int j = 0; j < nbVert; j++) {
			degree += poly.term(i).degree(A(0,j));
		}
#ifdef DEBUG
		cout << poly.term(i) << " Degree: " << degree;
#endif
		if(degree < maxDegree) {
			ex degreeUp = poly.term(i) * pow(basis, maxDegree - degree);
#ifdef DEBUG
			cout << "   --> New Term: " <<  degreeUp.expand();
#endif
			maxDegreePolynomial += degreeUp.expand();
		} else {
			maxDegreePolynomial += poly.term(i);
		}
#ifdef DEBUG
		cout << endl << "-----------------------------------------------" << endl;
#endif

	}
#ifdef DEBUG
	cout << endl << "Final Expansion: " << maxDegreePolynomial << endl << endl;
#endif
	return maxDegreePolynomial;
}



/*
 * Finds and prints the coefficients of the polynomial
 *
 *	maxDegreePolynomial: polynomial with monomials of degree n
 *	expandedBasis: basis powered n
 *	nbVert: number of vertices
 *	A: a_i matrix
 */
lst getCoefficients(ex &maxDegreePolynomial, polynomial &expandedBasis
		       , unsigned int nbVert, matrix &A)
{
	lst coefficients;

	cout << "-----------------------------------------------" << endl;
	cout << "Coefficients: " << endl << endl;

	for (size_t i = 0; i != expandedBasis.nbTerms(); ++i) {
		ex coeff = maxDegreePolynomial;
		ex vars = 1;
		for(unsigned int j = 0; j < nbVert; j++) {
			unsigned int deg = expandedBasis.term(i).degree(A(0,j));
			if(deg > 0) {
				coeff = coeff.coeff(A(0,j), deg);
				vars *= pow(A(0,j), deg);
			}
		}
		cout << "\t" << vars << "\t\t--> " << coeff / (expandedBasis.term(i)/vars);
		coefficients.append(coeff / (expandedBasis.term(i)/vars));
		cout << endl << endl;
	}
	return coefficients;
}


// replace the variables in the polynomial
ex replaceVariablesInPolynomial(ex &poly, const exvector &V, ex &variables)
{
	lst replace;

	for(unsigned int i = 0; i < V.size(); i++) {
#ifdef DEBUG
		cout << "Replacing: " << V(0,i) << " by " << variables[i] << endl;
#endif
		replace.append(V[i] == variables[i]);
	}
	ex polyRepl = poly.subs(replace);

	return(polyRepl.expand());
}


/* Converts int n to string */
string int2String(int n)
{
	char numeroVariable[DIGITS];
	snprintf(numeroVariable, DIGITS, "%d", n);
	string nroV(numeroVariable);

	return nroV;
}


/* 
 * Find the maximum multi-degree of the polinomial
 *
 *	polynomial: polynomial
 *	Vars: variables matrix
 *	nbVar: number of variables
 */
unsigned int findMaxDegree(ex polynomial, const exvector& Vars)
{
	unsigned int max = polynomial.degree(Vars[0]);

	for(unsigned int i = 0; i < Vars.size(); i++) {
		unsigned int degree = polynomial.degree(Vars[i]);
		if(max < degree) {
			max = degree;
		}
	}
	return max;
}
