#include <polylib/polylibgmp.h>
#include <stdlib.h>
#include <assert.h>
#include <util.h>

/* Return random value between 0 and max-1 inclusive
 */
int random_int(int max) {
    return (int) (((double)(max))*rand()/(RAND_MAX+1.0));
}

/* Inplace polarization
 */
void Polyhedron_Polarize(Polyhedron *P)
{
    unsigned NbRows = P->NbConstraints + P->NbRays;
    int i;
    Value **q;

    q = (Value **)malloc(NbRows * sizeof(Value *));
    assert(q);
    for (i = 0; i < P->NbRays; ++i)
	q[i] = P->Ray[i];
    for (; i < NbRows; ++i)
	q[i] = P->Constraint[i-P->NbRays];
    P->NbConstraints = NbRows - P->NbConstraints;
    P->NbRays = NbRows - P->NbRays;
    free(P->Constraint);
    P->Constraint = q;
    P->Ray = q + P->NbConstraints;
}

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
Polyhedron* supporting_cone(Polyhedron *P, int v)
{
    Matrix *M;
    Value tmp;
    int i, n, j;
    unsigned char *supporting = (unsigned char *)malloc(P->NbConstraints);
    unsigned dim = P->Dimension + 2;

    assert(v >=0 && v < P->NbRays);
    assert(value_pos_p(P->Ray[v][dim-1]));
    assert(supporting);

    value_init(tmp);
    for (i = 0, n = 0; i < P->NbConstraints; ++i) {
	Inner_Product(P->Constraint[i] + 1, P->Ray[v] + 1, dim - 1, &tmp);
	if ((supporting[i] = value_zero_p(tmp)))
	    ++n;
    }
    assert(n >= dim - 2);
    value_clear(tmp);
    M = Matrix_Alloc(n, dim);
    assert(M);
    for (i = 0, j = 0; i < P->NbConstraints; ++i)
	if (supporting[i]) {
	    value_set_si(M->p[j][dim-1], 0);
	    Vector_Copy(P->Constraint[i], M->p[j++], dim-1);
	}
    free(supporting);
    P = Constraints2Polyhedron(M, P->NbRays+1);
    assert(P);
    Matrix_Free(M);
    return P;
}

void value_lcm(Value i, Value j, Value* lcm)
{
    Value aux;
    value_init(aux);
    Gcd(i,j,lcm);
    value_multiply(aux,i,j);
    value_division(*lcm,aux,*lcm);
    value_clear(aux);
}

Polyhedron* supporting_cone_p(Polyhedron *P, Param_Vertices *v)
{
    Matrix *M;
    Value lcm, tmp;
    unsigned char *supporting = (unsigned char *)malloc(P->NbConstraints);
    unsigned dim = P->Dimension + 2;
    unsigned nparam = v->Vertex->NbColumns - 2;
    unsigned nvar = dim - nparam - 2;
    int i, n, j;
    Vector *row;

    assert(supporting);
    row = Vector_Alloc(nparam+1);
    assert(row);
    value_init(lcm);
    value_init(tmp);
    value_set_si(lcm, 1);
    for (i = 0, n = 0; i < P->NbConstraints; ++i) {
	Vector_Set(row->p, 0, nparam+1);
	for (j = 0 ; j < nvar; ++j) {
	    value_set_si(tmp, 1);
	    if (value_ne(lcm, v->Vertex->p[j][nparam+1])) {
		value_assign(tmp, lcm);
		value_lcm(lcm, v->Vertex->p[j][nparam+1], &lcm);
		value_division(tmp, lcm, tmp);
	    }
	    Vector_Combine(row->p, v->Vertex->p[j], row->p, tmp, P->Constraint[i][j+1], nparam+1);
	}
	value_set_si(tmp, 1);
	Vector_Combine(row->p, P->Constraint[i]+1+nvar, row->p, tmp, lcm, nparam+1);
	for (j = 0; j < nparam+1; ++j)
	    if (value_notzero_p(row->p[j]))
		break;
	if ((supporting[i] = (j == nparam + 1)))
	    ++n;
    }
    assert(n >= nvar);
    value_clear(tmp);
    value_clear(lcm);
    Vector_Free(row);
    M = Matrix_Alloc(n, nvar+2);
    assert(M);
    for (i = 0, j = 0; i < P->NbConstraints; ++i)
	if (supporting[i]) {
	    value_set_si(M->p[j][nvar+1], 0);
	    Vector_Copy(P->Constraint[i], M->p[j++], nvar+1);
	}
    free(supporting);
    P = Constraints2Polyhedron(M, P->NbRays+1);
    assert(P);
    Matrix_Free(M);
    return P;
}

Polyhedron* triangularize_cone(Polyhedron *P, unsigned NbMaxCons)
{
    const static int MAX_TRY=10;
    int i, j, r, n, t;
    Value tmp;
    unsigned dim = P->Dimension;
    Matrix *M = Matrix_Alloc(P->NbRays+1, dim+3);
    Matrix *M2;
    Polyhedron *L, *R, *T;
    assert(P->NbEq == 0);

    R = NULL;
    value_init(tmp);

    Vector_Set(M->p[0]+1, 0, dim+1);
    value_set_si(M->p[0][0], 1);
    value_set_si(M->p[0][dim+2], 1);
    Vector_Set(M->p[P->NbRays]+1, 0, dim+2);
    value_set_si(M->p[P->NbRays][0], 1);
    value_set_si(M->p[P->NbRays][dim+1], 1);

    for (i = 0, r = 1; i < P->NbRays; ++i) {
	if (value_notzero_p(P->Ray[i][dim+1]))
	    continue;
	Vector_Copy(P->Ray[i], M->p[r], dim+1);
	Inner_Product(M->p[r]+1, M->p[r]+1, dim, &tmp);
	value_assign(M->p[r][dim+1], tmp);
	value_set_si(M->p[r][dim+2], 0);
	++r;
    }

    L = Rays2Polyhedron(M, NbMaxCons);

    M2 = Matrix_Alloc(dim+1, dim+2);
    Vector_Set(M2->p[0]+1, 0, dim);
    value_set_si(M2->p[0][0], 1);
    value_set_si(M2->p[0][dim+1], 1);

    t = 1;
    if (0) {
try_again:
	/* Usually R should still be 0 */
	Domain_Free(R);
	Polyhedron_Free(L);
	for (r = 1; r < P->NbRays; ++r) {
	    value_set_si(M->p[r][dim+1], random_int((t+1)*dim)+1);
	}
	L = Rays2Polyhedron(M, NbMaxCons);
	++t;
    }
    assert(t <= MAX_TRY);

    R = NULL;
    n = 0;

    for (i = 0; i < L->NbConstraints; ++i) {
	if (value_negz_p(L->Constraint[i][dim+1]))
	    continue;
	if (value_notzero_p(L->Constraint[i][dim+2]))
	    continue;
	for (j = 1, r = 1; j < M->NbRows; ++j) {
	    Inner_Product(M->p[j]+1, L->Constraint[i]+1, dim+1, &tmp);
	    if (value_notzero_p(tmp))
		continue;
	    if (r > dim)
		goto try_again;
	    Vector_Copy(M->p[j]+1, M2->p[r]+1, dim);
	    value_set_si(M2->p[r][0], 1);
	    value_set_si(M2->p[r][dim+1], 0);
	    ++r;
	}
	assert(r == dim+1);
	T = Rays2Polyhedron(M2, P->NbConstraints+1);
	T->next = R;
	R = T;
	++n;
    }
    Matrix_Free(M2);

    Polyhedron_Free(L);
    value_clear(tmp);
    Matrix_Free(M);

    return R;
}

void check_triangulization(Polyhedron *P, Polyhedron *T)
{
    Polyhedron *C, *D, *E, *F, *G, *U;
    for (C = T; C; C = C->next) {
	if (C == T)
	    U = C;
	else 
	    U = DomainConvex(DomainUnion(U, C, 100), 100);
	for (D = C->next; D; D = D->next) {
	    F = C->next;
	    G = D->next;
	    C->next = NULL;
	    D->next = NULL;
	    E = DomainIntersection(C, D, 600);
	    assert(E->NbRays == 0 || E->NbEq >= 1);
	    Polyhedron_Free(E);
	    C->next = F;
	    D->next = G;
	}
    }
    assert(PolyhedronIncludes(U, P));
    assert(PolyhedronIncludes(P, U));
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
	value_set_si(*x, 0);
    else if (value_pos_p(a))
	value_assign(*x, e);
    else value_oppose(*x, e);
    if (value_zero_p(b))
	value_set_si(*y, 0);
    else {
	value_multiply(tmp, a, *x);
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
	value_assign(g, row->p[i]);
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
	assert(dim > 0);
	Vector_Gcd(p->Constraint[0]+1, dim+1, &g);
	Vector_AntiScale(p->Constraint[0]+1, p->Constraint[0]+1, g, dim+1);
	Vector_Gcd(p->Constraint[0]+1, dim, &g);
	if (value_notone_p(g) && value_notmone_p(g)) {
	    Polyhedron_Free(p);
	    p = Empty_Polyhedron(0);
	    break;
	}
	v = Vector_Alloc(dim);
	Vector_Copy(p->Constraint[0]+1, v->p, dim);
	m1 = unimodular_complete(v);
	m2 = Matrix_Alloc(dim, dim+1);
	for (i = 0; i < dim-1 ; ++i) {
	    Vector_Copy(m1->p[i+1], m2->p[i], dim);
	    value_set_si(m2->p[i][dim], 0);
	}
	Vector_Set(m2->p[dim-1], 0, dim);
	value_set_si(m2->p[dim-1][dim], 1);
	q = Polyhedron_Image(p, m2, p->NbConstraints+1+p->NbRays);
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

struct single {
    int	nr;
    int pos[2];
};

/*
 * The number of points in P is equal to factor time
 * the number of points in the polyhedron returned.
 * The return value is zero if no reduction can be found.
 */
Polyhedron* Polyhedron_Reduce(Polyhedron *P, Value* factor)
{
    int i, j, prev, nsingle, k, p;
    unsigned dim = P->Dimension;
    struct single *singles;
    int *bad;
    Value tmp, pos, neg;

    value_init(tmp);
    value_init(pos);
    value_init(neg);

    value_set_si(*factor, 1);

    singles = (struct single *)malloc(dim * sizeof(struct single));
    assert(singles);
    for (i = 0; i < dim; ++i)
	singles[i].nr = 0;
    bad = (int *)calloc(dim, sizeof(int));
    assert(bad);

    assert (P->NbEq == 0);

    nsingle = 0;
    for (i = 0; i < P->NbConstraints; ++i) {
	for (j = 0, prev = -1; j < dim; ++j) {
	    if (value_notzero_p(P->Constraint[i][j+1])) {
		if (prev == -1)
		    prev = j;
		else {
		    if (prev != -2)
			bad[prev] = 1;
		    bad[j] = 1;
		    prev = -2;
		}
	    }
	}
	if (prev >= 0)
	    singles[prev].pos[singles[prev].nr++] = i;
    }
    for (j = 0; j < dim; ++j) {
	if (bad[j])
	    singles[j].nr = 0;
	else if (singles[j].nr == 2)
	    ++nsingle;
    }
    if (nsingle) {
	Matrix *m = Matrix_Alloc((dim-nsingle)+1, dim+1);
	for (i = 0, j = 0; i < dim; ++i) {
	    if (singles[i].nr != 2)
		value_set_si(m->p[j++][i], 1);
	    else {
		for (k = 0; k <= 1; ++k) {
		    p = singles[i].pos[k];
		    value_oppose(tmp, P->Constraint[p][dim+1]);
		    if (value_pos_p(P->Constraint[p][i+1]))
			mpz_cdiv_q(pos, tmp, P->Constraint[p][i+1]);
		    else
			mpz_fdiv_q(neg, tmp, P->Constraint[p][i+1]);
		}
		value_substract(tmp, neg, pos);
		value_increment(tmp, tmp);
		value_multiply(*factor, *factor, tmp);
	    }
	}
	value_set_si(m->p[dim-nsingle][dim], 1);
	P = Polyhedron_Image(P, m, P->NbConstraints);
	Matrix_Free(m);
    } else
	P = NULL;
    free(singles);
    free(bad);

    value_clear(tmp);
    value_clear(pos);
    value_clear(neg);

    return P;
}

void manual_count(Polyhedron *P, Value* result)
{
    Polyhedron *U = Universe_Polyhedron(0);
    Enumeration *en = Polyhedron_Enumerate(P,U,1024);
    Value *v = compute_poly(en,NULL);
    value_assign(*result, *v);
    value_clear(*v);
    free(v);
    Enumeration_Free(en);
    Polyhedron_Free(U);
}

Bool isIdentity(Matrix *M)
{
    unsigned i, j;
    if (M->NbRows != M->NbColumns)
	return False;

    for (i = 0;i < M->NbRows; i ++)
	for (j = 0; j < M->NbColumns; j ++)
	    if (i == j) {
		if(value_notone_p(M->p[i][j]))
		    return False;
	    } else {
		if(value_notzero_p(M->p[i][j]))
		    return False;
	    }
    return True;
}

void Param_Polyhedron_Print(FILE* DST, Param_Polyhedron *PP, char **param_names)
{
  Param_Domain *P;
  Param_Vertices *V;

  for(P=PP->D;P;P=P->next) {
    
    /* prints current val. dom. */
    printf( "---------------------------------------\n" );
    printf( "Domain :\n");
    Print_Domain( stdout, P->Domain, param_names );
    
    /* scan the vertices */
    printf( "Vertices :\n");
    FORALL_PVertex_in_ParamPolyhedron(V,P,PP) {
	
      /* prints each vertex */
      Print_Vertex( stdout, V->Vertex, param_names );
      printf( "\n" );
    }
    END_FORALL_PVertex_in_ParamPolyhedron;
  }
}
