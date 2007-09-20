#include <assert.h>
#include <barvinok/util.h>
#include "remove_equalities.h"

static void transform(Polyhedron **P, Polyhedron **C, Matrix *CP, int free,
		      unsigned MaxRays)
{
    Polyhedron *Q = *P;
    Polyhedron *D = *C;
    Matrix *T;

    T = align_matrix(CP, Q->Dimension+1);
    *P = Polyhedron_Preimage(Q, T, MaxRays);
    if (free)
	Polyhedron_Free(Q);
    Matrix_Free(T);
    if (!D)
	return;
    *C = Polyhedron_Preimage(D, CP, MaxRays);
    if (free)
	Polyhedron_Free(D);
}

/* Remove all equalities in P and the context C (if not NULL).
 * Does not destroy P (or C).
 * Returns transformation on the parameters in the Matrix pointed to by CPP
 * (unless NULL) and transformation on the variables in the Matrix pointed to
 * by CVP (unless NULL).
 * Each of these transformation matrices maps the new parameters/variables
 * back to the original ones.
 */
int remove_all_equalities(Polyhedron **P, Polyhedron **C, Matrix **CPP, Matrix **CVP,
			  unsigned nparam, unsigned MaxRays)
{
    Matrix *CV = NULL;
    Matrix *CP = NULL;
    Polyhedron *Q = *P;
    Polyhedron *D = NULL;
    Polyhedron *R;
    int i;
    Matrix M;

    if (C) {
	D = *C;
	assert(D->Dimension == nparam);
    }

    if (Q->NbEq == 0 && (!D || D->NbEq == 0))
	return 0;

    if (D && D->NbEq) {
	Polyhedron_Matrix_View(D, &M, D->NbEq);
	CV = compress_variables(&M, 0);
	transform(&Q, &D, CV, Q != *P, MaxRays);
	nparam = CV->NbColumns-1;
    }

    /* compress_parms doesn't like equalities that only involve parameters */
    for (i = 0; i < Q->NbEq; ++i)
	if (First_Non_Zero(Q->Constraint[i]+1, Q->Dimension-nparam) == -1)
	    break;

    /* If we already compressed the parameters, then there should be
     * no such equalities left.
     */
    if (CV)
	assert(i >= Q->NbEq);

    if (i < Q->NbEq) {
	Matrix *M = Matrix_Alloc(Q->NbEq, 1+nparam+1);
	int n = 0;
	for (; i < Q->NbEq; ++i) {
	    if (First_Non_Zero(Q->Constraint[i]+1, Q->Dimension-nparam) == -1)
		Vector_Copy(Q->Constraint[i]+1+Q->Dimension-nparam,
			    M->p[n++]+1, nparam+1);
	}
	M->NbRows = n;
	CV = compress_variables(M, 0);
	Matrix_Free(M);
	transform(&Q, &D, CV, Q != *P, MaxRays);
	nparam = CV->NbColumns-1;
    }

    if (emptyQ2(Q)) {
	if (CV)
	    Matrix_Free(CV);
	*P = Q;
	if (C)
	    *C = D;
	return 0;
    }

    if (Q->NbEq == 0) {
	CP = CV;
	CV = NULL;
    } else {
	Polyhedron_Matrix_View(Q, &M, Q->NbEq);
	CP = compress_parms(&M, nparam);

	if (isIdentity(CP)) {
	    Matrix_Free(CP);
	    CP = CV;
	    CV = NULL;
	} else {
	    transform(&Q, &D, CP, Q != *P, MaxRays);
	    if (CV) {
		Matrix *T = CP;
		CP = Matrix_Alloc(CV->NbRows, T->NbColumns);
		Matrix_Product(CV, T, CP);
		Matrix_Free(T);
		Matrix_Free(CV);
		CV = NULL;
	    }
	    nparam = CP->NbColumns-1;
	}

	Polyhedron_Matrix_View(Q, &M, Q->NbEq);
	CV = compress_variables(&M, nparam);
	if (!CV) {
	    if (Q != *P)
		Polyhedron_Free(Q);
	    Q = NULL;
	} else if (isIdentity(CV)) {
	    Matrix_Free(CV);
	    CV = NULL;
	} else {
	    R = Polyhedron_Preimage(Q, CV, MaxRays);
	    if (Q != *P)
		Polyhedron_Free(Q);
	    Q = R;
	}
    }

    *P = Q;
    if (C)
	*C = D;
    if (CPP)
	*CPP = CP;
    else if (CP)
	Matrix_Free(CP);
    if (CVP)
	*CVP = CV;
    else if (CV)
	Matrix_Free(CV);
    return 1;
}
