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

    L = Rays2Polyhedron(M, NbMaxCons);

    if (L->NbEq != 0) {
	Polyhedron_Free(L);
	value_set_si(tmp, 2);
	Vector_Scale(M->p[1]+1, M->p[1]+1, tmp, dim);
	Inner_Product(M->p[1]+1, M->p[1]+1, dim, &tmp);
	value_assign(M->p[1][dim+1], tmp);
	L = Rays2Polyhedron(M, NbMaxCons);
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
	T = Rays2Polyhedron(M2, P->NbConstraints);
	T->next = R;
	R = T;
    }
    Matrix_Free(M2);

    Polyhedron_Free(L);
    value_clear(tmp);
    Matrix_Free(M);

    return R;
}

void Euclid(Value a, Value b, Value *x, Value *y, Value *g)
{
    Value c, d, e, f, tmp;

    value_init(c);
    value_init(d);
    value_init(e);
    value_init(f);
    value_init(tmp);
    value_absolute(c, a);
    value_absolute(d, b);
    value_set_si(e, 1);
    value_set_si(f, 0);
    while(value_pos_p(d)) {
	value_division(tmp, c, d);
	value_multiply(tmp, tmp, f);
	value_substract(e, e, tmp);
	value_division(tmp, c, d);
	value_multiply(tmp, tmp, d);
	value_substract(c, c, tmp);
	value_swap(c, d);
	value_swap(e, f);
    }
    value_assign(*g, c);
    if (value_zero_p(a))
	value_assign(*x, 0);
    else if (value_pos_p(a))
	value_assign(*x, e);
    else value_oppose(*x, e);
    if (value_zero_p(b))
	value_assign(*y, 0);
    else {
	value_multiply(tmp, a, e);
	value_substract(tmp, c, tmp);
	value_division(*y, tmp, b);
    }
    value_clear(c);
    value_clear(d);
    value_clear(e);
    value_clear(f);
    value_clear(tmp);
}

Matrix * unimodular_complete(Vector *row) 
{
    Value g, b, c, old, tmp;
    Matrix *m;
    unsigned i, j;

    value_init(b);
    value_init(c);
    value_init(g);
    value_init(old);
    value_init(tmp);
    m = Matrix_Alloc(row->Size, row->Size);
    for (j = 0; j < row->Size; ++j) {
	value_assign(m->p[0][j], row->p[j]);
    }
    value_assign(g, row->p[0]);
    for (i = 1; value_zero_p(g) && i < row->Size; ++i) {
	for (j = 0; j < row->Size; ++j) {
	    if (j == i-1)
		value_set_si(m->p[i][j], 1);
	    else
		value_set_si(m->p[i][j], 0);
	}
    }
    for (; i < row->Size; ++i) {
	value_assign(old, g);
	Euclid(old, row->p[i], &c, &b, &g);
	value_oppose(b, b);
	for (j = 0; j < row->Size; ++j) {
	    if (j < i) {
		value_multiply(tmp, row->p[j], b);
		value_division(m->p[i][j], tmp, old);
	    } else if (j == i)
		value_assign(m->p[i][j], c);
	    else
		value_set_si(m->p[i][j], 0);
	}
    }
    value_clear(b);
    value_clear(c);
    value_clear(g);
    value_clear(old);
    value_clear(tmp);
    return m;
}

/*
 * Returns a full-dimensional polyhedron with the same number
 * of integer points as P
 */
Polyhedron *remove_equalities(Polyhedron *P)
{
    Value g;
    Vector *v;
    Polyhedron *p = Polyhedron_Copy(P), *q;
    unsigned dim = p->Dimension;
    Matrix *m1, *m2;
    int i;

    value_init(g);
    while (p->NbEq > 0) {
	v = Vector_Alloc(dim);
	Vector_Gcd(p->Constraint[0]+1, dim, &g);
	Vector_AntiScale(p->Constraint[0]+1, v->p, g, dim);
	m1 = unimodular_complete(v);
	m2 = Matrix_Alloc(dim, dim+1);
	for (i = 0; i < dim-1 ; ++i) {
	    Vector_Copy(m1->p[i+1], m2->p[i], dim);
	    value_set_si(m2->p[i][dim], 0);
	}
	Vector_Set(m2->p[dim-1], 0, dim);
	value_set_si(m2->p[dim-1][dim], 1);
	q = Polyhedron_Image(p, m2, p->NbConstraints);
	Vector_Free(v);
	Matrix_Free(m1);
	Matrix_Free(m2);
	Polyhedron_Free(p);
	p = q;
	--dim;
    }
    value_clear(g);
    return p;
}

void manual_count(Polyhedron *P, Value* result)
{
    Polyhedron *U = Universe_Polyhedron(0);
    Enumeration *ee, *en = Polyhedron_Enumerate(P,U,1024);
    Value *v = compute_poly(en,NULL);
    value_assign(*result, *v);
    value_clear(*v);
    free(v);
    Enumeration_Free(en);
    Polyhedron_Free(U);
}
