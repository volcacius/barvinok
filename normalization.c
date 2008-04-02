#include <assert.h>
#include "normalization.h"

static int is_unit_row(Value *row, int pos, int len)
{
    if (!value_one_p(row[pos]) && !value_mone_p(row[pos]))
	return 0;
    return First_Non_Zero(row+pos+1, len-(pos+1)) == -1;
}

/* Transform the constraints C to "standard form".
 * In particular, if C is described by
 *		A x + b(p) >= 0
 * then this function returns a matrix H = A U, A = H Q, such
 * that D x' = D Q x >= -b(p), with D a diagonal matrix with
 * positive entries.  The calling function can then construct
 * the standard form H' x' - I s + b(p) = 0, with H' the rows of H
 * that are not positive multiples of unit vectors
 * (since those correspond to D x' >= -b(p)).
 * The number of rows in H' is returned in *rows_p.
 * Optionally returns the matrix that maps the new variables
 * back to the old variables x = U x'.
 * Note that the rows of H (and C) may be reordered by this function.
 */
Matrix *standard_constraints(Matrix *C, unsigned nparam, int *rows_p,
			     Matrix **T)
{
    unsigned dim = C->NbColumns - 2;
    unsigned nvar = dim - nparam;
    int i, j, d;
    int rows;
    Matrix *M;
    Matrix *H, *U, *Q;

    for (i = 0; i < C->NbRows; ++i)
	assert(value_one_p(C->p[i][0]));

    /* move constraints only involving parameters down
     * and move unit vectors (if there are any) to the right place.
     */
    for (d = 0, j = C->NbRows; d < j; ++d) {
	int index;
	index = First_Non_Zero(C->p[d]+1, nvar);
	if (index != -1) {
	    if (index != d &&
		is_unit_row(C->p[d]+1, index, nvar)) {
		Vector_Exchange(C->p[d], C->p[index], dim+2);
		if (index > d)
		    --d;
	    }
	    continue;
	}
	while (d < --j && First_Non_Zero(C->p[j]+1, nvar) == -1)
	    ;
	if (d >= j)
	    break;
	Vector_Exchange(C->p[d], C->p[j], dim+2);
    }
    M = Matrix_Alloc(d, nvar);
    for (j = 0; j < d; ++j)
	Vector_Copy(C->p[j]+1, M->p[j], nvar);

    neg_left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    Matrix_Free(Q);
    if (T)
	*T = U;
    else
	Matrix_Free(U);

    /* Rearrange rows such that top of H is lower diagonal and
     * compute the number of non (multiple of) unit-vector rows.
     */
    rows = H->NbRows-nvar;
    for (i = 0; i < H->NbColumns; ++i) {
	for (j = i; j < H->NbRows; ++j)
	    if (value_notzero_p(H->p[j][i]))
		break;
	if (j != i) {
	    Vector_Exchange(C->p[i], C->p[j], dim+2);
	    Vector_Exchange(H->p[i], H->p[j], H->NbColumns);
	}
	if (First_Non_Zero(H->p[i], i) != -1)
	    rows++;
    }
    *rows_p = rows;

    return H;
}
