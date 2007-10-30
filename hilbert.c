#define Vector ZSolveVector
#define Matrix ZSolveMatrix
#include "zsolve/libzsolve.h"
#undef Vector
#undef Matrix
#include "hilbert.h"

static ZSolveMatrix Matrix2zsolve(Matrix *M)
{
    int i, j;
    ZSolveMatrix zmatrix;

    zmatrix = createMatrix(M->NbColumns-2, M->NbRows);
    for (i = 0; i < M->NbRows; ++i)
	for (j = 0; j < M->NbColumns-2; ++j) {
	    assert(mpz_cmp_si(M->p[i][1+j], -MAXINT) > 0);
	    assert(mpz_cmp_si(M->p[i][1+j], MAXINT) < 0);
	    zmatrix->Data[i*zmatrix->Width+j] = mpz_get_si(M->p[i][1+j]);
	}

    return zmatrix;
}

static Matrix *VectorArray2Matrix(VectorArray array, unsigned cols)
{
    int i, j;
    Matrix *M = Matrix_Alloc(array->Size, cols+1);

    for (i = 0; i < array->Size; ++i) {
	for (j = 0; j < cols; ++j)
	    value_set_si(M->p[i][j], array->Data[i][j]);
	value_set_si(M->p[i][cols], 1);
    }
    return M;
}

static void Polyhedron_Remove_Positivity_Constraint(Polyhedron *P)
{
    int i;

    for (i = 0; i < P->NbConstraints; ++i) {
	if (First_Non_Zero(P->Constraint[i]+1, P->Dimension) != -1)
	    continue;
	if (i < P->NbConstraints-1)
	    Vector_Exchange(P->Constraint[i],
			    P->Constraint[P->NbConstraints-1],
			    P->Dimension+2);
	P->NbConstraints--;
	--i;
    }
}

static Matrix *Polyhedron2standard_form(Polyhedron *P, Matrix **T)
{
    int i;
    unsigned dim = P->Dimension;
    Matrix *M2 = Matrix_Alloc(P->NbConstraints+1, dim+1);
    Matrix *H, *Q;

    Polyhedron_Remove_Positivity_Constraint(P);
    for (i = 0; i < P->NbConstraints; ++i) {
	assert(value_zero_p(P->Constraint[i][1+dim]));
	Vector_Copy(P->Constraint[i]+1, M2->p[i], dim);
    }
    value_set_si(M2->p[P->NbConstraints][dim], 1);
    neg_left_hermite(M2, &H, &Q, T);
    Matrix_Free(Q);
    Matrix_Free(M2);

    M2 = Matrix_Alloc(P->NbConstraints, 2+dim+(P->NbConstraints-P->NbEq));

    for (i = 0; i < P->NbEq; ++i)
	Vector_Copy(H->p[i], M2->p[i]+1, dim);
    for (i = P->NbEq; i < P->NbConstraints; ++i) {
	Vector_Copy(H->p[i], M2->p[i]+1, dim);
	value_set_si(M2->p[i][1+dim+i-P->NbEq], -1);
    }
    Matrix_Free(H);
    return M2;
}

/* Assumes C is a linear cone (i.e. with apex zero).
 * All equalities are removed first to speed up the computation
 * in zsolve.
 */
Matrix *Cone_Hilbert_Basis(Polyhedron *C, unsigned MaxRays)
{
    unsigned dim;
    int i;
    Matrix *M2, *M3, *T;
    Matrix *CV = NULL;
    LinearSystem initialsystem;
    ZSolveMatrix matrix;
    ZSolveVector rhs;
    ZSolveContext ctx;

    remove_all_equalities(&C, NULL, NULL, &CV, 0, MaxRays);
    dim = C->Dimension;

    for (i = 0; i < C->NbConstraints; ++i)
	assert(value_zero_p(C->Constraint[i][1+dim]) ||
	       First_Non_Zero(C->Constraint[i]+1, dim) == -1);

    M2 = Polyhedron2standard_form(C, &T);
    matrix = Matrix2zsolve(M2);
    Matrix_Free(M2);

    rhs = createVector(matrix->Height);
    for (i = 0; i < matrix->Height; i++)
	rhs[i] = 0;

    initialsystem = createLinearSystem();
    setLinearSystemMatrix(initialsystem, matrix);
    deleteMatrix(matrix);

    setLinearSystemRHS(initialsystem, rhs);
    deleteVector(rhs);

    setLinearSystemLimit(initialsystem, -1, 0, MAXINT, 0);
    setLinearSystemEquationType(initialsystem, -1, EQUATION_EQUAL, 0);

    ctx = createZSolveContextFromSystem(initialsystem, NULL, 0, 0, NULL, NULL);
    zsolveSystem(ctx, 0);

    M2 = VectorArray2Matrix(ctx->Homs, C->Dimension);
    deleteZSolveContext(ctx, 1);
    Matrix_Transposition(T);
    M3 = Matrix_Alloc(M2->NbRows, M2->NbColumns);
    Matrix_Product(M2, T, M3);
    Matrix_Free(M2);
    Matrix_Free(T);

    if (CV) {
	Matrix *T, *M;
	T = Transpose(CV);
	M = Matrix_Alloc(M3->NbRows, T->NbColumns);
	Matrix_Product(M3, T, M);
	Matrix_Free(M3);
	Matrix_Free(CV);
	Matrix_Free(T);
	Polyhedron_Free(C);
	M3 = M;
    }

    return M3;
}
