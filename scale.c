#include <barvinok/util.h>
#include "scale.h"

/* If a vertex is described by A x + B p + c = 0, then
 * M = [A B] and we want to compute a linear transformation L such
 * that H L = A and H \Z contains both A \Z and B \Z.
 * We compute
 *             [ A B ] = [ H 0 ] [ U_11  U_12 ]
 *                               [ U_21  U_22 ]
 *
 * U_11 is the required linear transformation.
 * Note that this also works if M has more rows than there are variables,
 * i.e., if some rows in M are linear combinations of other rows.
 * These extra rows only affect and H and not U.
 */
static Lattice *extract_lattice(Matrix *M, unsigned nvar)
{
    int row;
    Matrix *H, *Q, *U, *Li;
    Lattice *L;
    int ok;

    left_hermite(M, &H, &Q, &U);
    Matrix_Free(U);

    Li = Matrix_Alloc(nvar+1, nvar+1);
    L = Matrix_Alloc(nvar+1, nvar+1);
    value_set_si(Li->p[nvar][nvar], 1);

    for (row = 0; row < nvar; ++row)
	Vector_Copy(Q->p[row], Li->p[row], nvar);
    Matrix_Free(H);
    Matrix_Free(Q);

    ok = Matrix_Inverse(Li, L);
    assert(ok);
    Matrix_Free(Li);

    return L;
}

/* Returns the smallest (wrt inclusion) lattice that contains both X and Y */
static Lattice *LatticeJoin(Lattice *X, Lattice *Y)
{
    int i;
    int dim = X->NbRows-1;
    Value lcm;
    Value tmp;
    Lattice *L;
    Matrix *M, *H, *U, *Q;

    assert(X->NbColumns-1 == dim);
    assert(Y->NbRows-1 == dim);
    assert(Y->NbColumns-1 == dim);

    value_init(lcm);
    value_init(tmp);

    M = Matrix_Alloc(dim, 2*dim);
    value_lcm(X->p[dim][dim], Y->p[dim][dim], &lcm);

    value_division(tmp, lcm, X->p[dim][dim]);
    for (i = 0; i < dim; ++i)
	Vector_Scale(X->p[i], M->p[i], tmp, dim);
    value_division(tmp, lcm, Y->p[dim][dim]);
    for (i = 0; i < dim; ++i)
	Vector_Scale(Y->p[i], M->p[i]+dim, tmp, dim);

    left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    Matrix_Free(Q);
    Matrix_Free(U);

    L = Matrix_Alloc(dim+1, dim+1);
    value_assign(L->p[dim][dim], lcm);
    for (i = 0; i < dim; ++i)
	Vector_Copy(H->p[i], L->p[i], dim);
    Matrix_Free(H);

    value_clear(tmp);
    value_clear(lcm);
    return L;
}

static void Param_Vertex_Common_Denominator(Param_Vertices *V)
{
    unsigned dim;
    Value lcm;
    int i;

    assert(V->Vertex->NbRows > 0);
    dim = V->Vertex->NbColumns-2;

    value_init(lcm);

    value_assign(lcm, V->Vertex->p[0][dim+1]);
    for (i = 1; i < V->Vertex->NbRows; ++i)
	value_lcm(V->Vertex->p[i][dim+1], lcm, &lcm);

    for (i = 0; i < V->Vertex->NbRows; ++i) {
	if (value_eq(V->Vertex->p[i][dim+1], lcm))
	    continue;
	value_division(V->Vertex->p[i][dim+1], lcm, V->Vertex->p[i][dim+1]);
	Vector_Scale(V->Vertex->p[i], V->Vertex->p[i],
		     V->Vertex->p[i][dim+1], dim+1);
	value_assign(V->Vertex->p[i][dim+1], lcm);
    }

    value_clear(lcm);
}

static void Param_Vertex_Image(Param_Vertices *V, Matrix *T)
{
    unsigned nvar  = V->Vertex->NbRows;
    unsigned nparam = V->Vertex->NbColumns - 2;
    Matrix *Vertex;
    int i;

    Param_Vertex_Common_Denominator(V);
    Vertex = Matrix_Alloc(V->Vertex->NbRows, V->Vertex->NbColumns);
    Matrix_Product(T, V->Vertex, Vertex);
    for (i = 0; i < nvar; ++i) {
	value_assign(Vertex->p[i][nparam+1], V->Vertex->p[i][nparam+1]);
	Vector_Normalize(Vertex->p[i], nparam+2);
    }
    Matrix_Free(V->Vertex);
    V->Vertex = Vertex;
}

/* Scales the parametric polyhedron with constraints *P and vertices PP
 * such that the number of integer points can be represented by a polynomial.
 * Both *P and P->Vertex are adapted according to the scaling.
 * The scaling factor is returned in *det.
 * The enumerator of the scaled parametric polyhedron should be divided
 * by this number to obtain an approximation of the enumerator of the
 * original parametric polyhedron.
 *
 * The algorithm is described in "Approximating Ehrhart Polynomials using
 * affine transformations" by B. Meister.
 */
void Param_Polyhedron_Scale_Integer_Slow(Param_Polyhedron *PP, Polyhedron **P,
					 Value *det, unsigned MaxRays)
{
    Param_Vertices *V;
    unsigned dim = (*P)->Dimension;
    unsigned nparam;
    unsigned nvar;
    Lattice *L = NULL, *Li;
    Matrix *T;
    Matrix *expansion;
    int i;
    int ok;

    if (!PP->nbV)
	return;

    nparam = PP->V->Vertex->NbColumns - 2;
    nvar = dim - nparam;

    for (V = PP->V; V; V = V->next) {
	Lattice *L2;
	Matrix *M;
	int i, j, n;
	unsigned char *supporting;

	supporting = supporting_constraints(*P, V, &n);
	M = Matrix_Alloc(n, (*P)->Dimension);
	for (i = 0, j = 0; i < (*P)->NbConstraints; ++i)
	    if (supporting[i])
		Vector_Copy((*P)->Constraint[i]+1, M->p[j++], (*P)->Dimension);
	free(supporting);
	L2 = extract_lattice(M, nvar);
	Matrix_Free(M);

	if (!L)
	    L = L2;
	else {
	    Lattice *L3 = LatticeJoin(L, L2);
	    Matrix_Free(L);
	    Matrix_Free(L2);
	    L = L3;
	}
    }

    /* apply the variable expansion to the polyhedron (constraints) */
    expansion = Matrix_Alloc(nvar + nparam + 1,  nvar + nparam + 1);
    for (i = 0; i < nvar; ++i)
	Vector_Copy(L->p[i], expansion->p[i], nvar);
    for (i = nvar; i < nvar+nparam+1; ++i)
	value_assign(expansion->p[i][i], L->p[nvar][nvar]);

    *P = Polyhedron_Preimage(*P, expansion, MaxRays);
    Matrix_Free(expansion);

    /* apply the variable expansion to the parametric vertices */
    Li = Matrix_Alloc(nvar+1, nvar+1);
    ok = Matrix_Inverse(L, Li);
    assert(ok);
    Matrix_Free(L);
    assert(value_one_p(Li->p[nvar][nvar]));
    T = Matrix_Alloc(nvar, nvar);
    value_set_si(*det, 1);
    for (i = 0; i < nvar; ++i) {
	value_multiply(*det, *det, Li->p[i][i]);
	Vector_Copy(Li->p[i], T->p[i], nvar);
    }
    Matrix_Free(Li);
    for (V = PP->V; V; V = V->next)
	Param_Vertex_Image(V, T);
    Matrix_Free(T);
}
