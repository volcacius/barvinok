#include <polylib/polylibgmp.h>
#include <assert.h>

/*
 * Rather general polar
 * We can optimize it significantly if we assume that
 * P includes zero
 *
 * Also, we calculate the polar as defined in Schrijver
 * The opposite should probably work as well and would
 * eliminate the need for multiplying by -1
 */
Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays)
{
    int i;
    Value mone;
    unsigned dim = P->Dimension + 2;
    Matrix *M = Matrix_Alloc(P->NbRays, dim);

    assert(M);
    value_init(mone);
    value_set_si(mone, -1);
    for (i = 0; i < P->NbRays; ++i) {
	Vector_Scale(P->Ray[i], M->p[i], mone, dim);
	value_multiply(M->p[i][0], M->p[i][0], mone);
	value_multiply(M->p[i][dim-1], M->p[i][dim-1], mone);
    }
    P = Constraints2Polyhedron(M, NbMaxRays);
    assert(P);
    Matrix_Free(M);
    value_clear(mone);
    return P;
}

/*
 * Returns the supporting cone of P at the vertex with index v
 */
Polyhedron* supporting_cone(Polyhedron *P, int v, unsigned NbMaxRays)
{
    Matrix *M;
    Value tmp;
    int i, n, j;
    unsigned char *supporting = (unsigned char *)malloc(P->NbRays);
    unsigned dim = P->Dimension + 2;

    assert(v >=0 && v < P->NbRays);
    assert(value_pos_p(P->Ray[v][dim-1]));
    assert(supporting);

    value_init(tmp);
    for (i = 0, n = 0; i < P->NbRays; ++i) {
	Inner_Product(P->Constraint[i] + 1, P->Ray[v] + 1, dim - 1, &tmp);
	if ((supporting[i] = value_zero_p(tmp)))
	    ++n;
    }
    value_clear(tmp);
    M = Matrix_Alloc(n, dim);
    assert(M);
    for (i = 0, j = 0; i < P->NbRays; ++i)
	if (supporting[i])
	    Vector_Copy(P->Constraint[i], M->p[j++], dim);
    free(supporting);
    P = Constraints2Polyhedron(M, NbMaxRays);
    assert(P);
    Matrix_Free(M);
    M = Matrix_Alloc(P->NbRays, dim);
    assert(M);
    n = P->NbConstraints;
    Vector_Copy(P->Ray[0], M->p[0], P->NbRays * dim);
    Polyhedron_Free(P);
    for (i = 0; i < M->NbRows; ++i)
	if (value_notzero_p(M->p[i][dim-1])) {
	    Vector_Set(M->p[i]+1, 0, dim - 2);
	    break;
	}
    P = Rays2Polyhedron(M, n+1);
    assert(P);
    Matrix_Free(M);
    return P;
}
