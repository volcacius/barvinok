#define Vector ZSolveVector
#define Matrix ZSolveMatrix
#include "zsolve/libzsolve.h"
#undef Vector
#undef Matrix
#include "hilbert.h"
#include "topcom.h"

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

/* Return
 *             T 0
 *             0 1
 */
static Matrix *homogenize(Matrix *T)
{
    int i;
    Matrix *H = Matrix_Alloc(T->NbRows+1, T->NbColumns+1);

    for (i = 0; i < T->NbRows; ++i)
	Vector_Copy(T->p[i], H->p[i], T->NbColumns);
    value_set_si(H->p[T->NbRows][T->NbColumns], 1);
    return H;
}

static Matrix *Polyhedron2standard_form(Polyhedron *P, Matrix **T)
{
    int i, j;
    int rows;
    unsigned dim = P->Dimension;
    Matrix *M2;
    Matrix *H, *U;

    assert(P->NbEq == 0);
    Polyhedron_Remove_Positivity_Constraint(P);
    for (i = 0; i < P->NbConstraints; ++i)
	assert(value_zero_p(P->Constraint[i][1+dim]));

    H = standard_constraints(P, 0, &rows, &U);
    *T = homogenize(U);
    Matrix_Free(U);

    M2 = Matrix_Alloc(rows, 2+dim+rows);

    for (i = dim; i < H->NbRows; ++i) {
	Vector_Copy(H->p[i], M2->p[i-dim]+1, dim);
	value_set_si(M2->p[i-dim][1+i], -1);
    }
    for (i = 0, j = H->NbRows-dim; i < dim; ++i) {
	if (First_Non_Zero(H->p[i], i) == -1)
	    continue;
	Vector_Oppose(H->p[i], M2->p[j]+1, dim);
	value_set_si(M2->p[j][1+j+dim], 1);
	++j;
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
