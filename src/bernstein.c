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
			Matrix_Set(i, j, (int) M[i*columns+j], retval);
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
