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

Polyhedron* triangularize_cone(Polyhedron *P, unsigned NbMaxCons)
{
    int i, j, r;
    Value tmp;
    unsigned dim = P->Dimension;
    Matrix *M = Matrix_Alloc(P->NbRays, dim+3);
    Matrix *M2;
    Polyhedron *L, *R, *T;

    R = NULL;
    value_init(tmp);

    Vector_Set(M->p[0]+1, 0, dim+1);
    value_set_si(M->p[0][0], 1);
    value_set_si(M->p[0][dim+2], 1);

    for (i = 0, r = 1; i < P->NbRays; ++i) {
	if (value_notzero_p(P->Ray[i][dim+1]))
	    continue;
	Vector_Copy(P->Ray[i], M->p[r], dim+1);
	Inner_Product(M->p[r]+1, M->p[r]+1, dim, &tmp);
	value_assign(M->p[r][dim+1], tmp);
	value_set_si(M->p[r][dim+2], 1);
	++r;
    }

    Matrix_Print(stdout, P_VALUE_FMT, M);
    L = Rays2Polyhedron(M, NbMaxCons);
    Polyhedron_Print(stdout, P_VALUE_FMT, L);

    if (L->NbEq != 0) {
	Polyhedron_Free(L);
	value_set_si(tmp, 2);
	Vector_Scale(M->p[1]+1, M->p[1]+1, tmp, dim);
	Inner_Product(M->p[1]+1, M->p[1]+1, dim, &tmp);
	value_assign(M->p[1][dim+1], tmp);
	Matrix_Print(stdout, P_VALUE_FMT, M);
	L = Rays2Polyhedron(M, NbMaxCons);
	Polyhedron_Print(stdout, P_VALUE_FMT, L);
	assert(L->NbEq == 0);
    }

    M2 = Matrix_Alloc(dim+1, dim+2);
    Vector_Set(M2->p[0]+1, 0, dim);
    value_set_si(M2->p[0][0], 1);
    value_set_si(M2->p[0][dim+1], 1);
    for (i = 0; i < L->NbConstraints; ++i) {
	assert(value_notzero_p(L->Constraint[i][dim+1]));
	if (value_neg_p(L->Constraint[i][dim+1]))
	    continue;
	if (value_notzero_p(L->Constraint[i][dim+2]))
	    continue;
	for (j = 1, r = 1; j < M->NbRows; ++j) {
	    Inner_Product(M->p[j]+1, L->Constraint[i]+1, dim+1, &tmp);
	    if (value_notzero_p(tmp))
		continue;
	    Vector_Copy(M->p[j]+1, M2->p[r]+1, dim);
	    value_set_si(M2->p[r][0], 1);
	    value_set_si(M2->p[r][dim+1], 0);
	    ++r;
	}
	assert(r == dim+1);
	Matrix_Print(stdout, P_VALUE_FMT, M2);
	T = Rays2Polyhedron(M2, P->NbConstraints);
	Polyhedron_Print(stdout, P_VALUE_FMT, T);
	T->next = R;
	R = T;
    }
    Matrix_Free(M2);

    Polyhedron_Free(L);
    value_clear(tmp);
    Matrix_Free(M);

    return R;
}
