/* 
 *	Bernstein Expansion
 *
 *	- polylib functions
 */


#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include <polylib/polylibgmp.h>

#include "bernstein.h"

// global to be used on the coefficient max
extern Polyhedron *VD = NULL;	// current validity domain

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



/* Prints the domain corresponding to Q */
void printDomain(Param_Domain *Q, Polyhedron *B, char **param_name)
{

	printf("\nDomain: \n");

	if(VD != NULL) {
		Domain_Free(VD);
	}
	VD = DomainSimplify(Q->Domain, B, MAXRAYS);
	Print_Domain(stdout, VD, param_name);
}


/* Converts a Matrix from the polylib format to long long* format */
long long *matrix2longlong(Matrix *M)
{
 	long long *matrix = calloc(sizeof(long long), M->NbRows * M->NbColumns);
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





/* Sets a value on the Matrix A, A(r, c) = v */
void Matrix_Set(int r, int c, int v, Matrix * A)
{
        Value val;
        value_init(val);
        value_set_si(val,v);
        value_assign(A->p[r][c] , val);
}

/* Converts a *longlong matrix to polylib format */
Matrix *longlong2polylib(long long *M, unsigned int rows, unsigned int columns)
{
	Matrix *retval;
	unsigned int i, j;

	retval = Matrix_Alloc(rows, columns);
	for(i = 0; i < rows; i++) {
		for(j = 0; j < columns; j++) {
			Matrix_Set(i, j, (int) M[i*rows+j], retval);
		}
	}
	return retval;
}


/* Check  if a given set of constraints (M matrix) holds or not in the validity domain */
unsigned  checkConstraint(long long *M, unsigned int rows, unsigned int columns)
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
	newB = DomainIntersection(VD, mC, MAXRAYS);

#ifdef DEBUG
	printf("Intersection betwenn Max/Min/Pos/Neg Constraints and \nOriginal Validity Domain: \n");
	Print_Domain(stdout, newB, param_name);
#endif

	if(emptyQ(newB)) {
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
