#include <polylib/polylibgmp.h>
#include <stdlib.h>
#include <assert.h>
#include "config.h"

#ifndef HAVE_ENUMERATE4
#define Polyhedron_Enumerate(a,b,c,d) Polyhedron_Enumerate(a,b,c)
#endif

void manual_count(Polyhedron *P, Value* result)
{
    Polyhedron *U = Universe_Polyhedron(0);
    Enumeration *en = Polyhedron_Enumerate(P,U,1024,NULL);
    Value *v = compute_poly(en,NULL);
    value_assign(*result, *v);
    value_clear(*v);
    free(v);
    Enumeration_Free(en);
    Polyhedron_Free(U);
}

#include "ev_operations.h"
#include <util.h>
#include <barvinok.h>

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
    value_multiply(aux,i,j);
    Gcd(i,j,lcm);
    value_division(*lcm,aux,*lcm);
    value_clear(aux);
}

Polyhedron* supporting_cone_p(Polyhedron *P, Param_Vertices *v)
{
    Matrix *M;
    Value lcm, tmp, tmp2;
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
    value_init(tmp2);
    value_set_si(lcm, 1);
    for (i = 0, n = 0; i < P->NbConstraints; ++i) {
	Vector_Set(row->p, 0, nparam+1);
	for (j = 0 ; j < nvar; ++j) {
	    value_set_si(tmp, 1);
	    value_assign(tmp2,  P->Constraint[i][j+1]);
	    if (value_ne(lcm, v->Vertex->p[j][nparam+1])) {
		value_assign(tmp, lcm);
		value_lcm(lcm, v->Vertex->p[j][nparam+1], &lcm);
		value_division(tmp, lcm, tmp);
		value_multiply(tmp2, tmp2, lcm);
		value_division(tmp2, tmp2, v->Vertex->p[j][nparam+1]);
	    }
	    Vector_Combine(row->p, v->Vertex->p[j], row->p, 
			   tmp, tmp2, nparam+1);
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
    value_clear(tmp2);
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
    Matrix *M2, *M3;
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

    M3 = Matrix_Copy(M);
    L = Rays2Polyhedron(M3, NbMaxCons);
    Matrix_Free(M3);

    M2 = Matrix_Alloc(dim+1, dim+2);

    t = 1;
    if (0) {
try_again:
	/* Usually R should still be 0 */
	Domain_Free(R);
	Polyhedron_Free(L);
	for (r = 1; r < P->NbRays; ++r) {
	    value_set_si(M->p[r][dim+1], random_int((t+1)*dim)+1);
	}
	M3 = Matrix_Copy(M);
	L = Rays2Polyhedron(M3, NbMaxCons);
	Matrix_Free(M3);
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
	Vector_Set(M2->p[0]+1, 0, dim);
	value_set_si(M2->p[0][0], 1);
	value_set_si(M2->p[0][dim+1], 1);
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

/*
 * Returns a full-dimensional polyhedron with the same number
 * of integer points as P
 * nvar specifies the number of variables
 * The remaining dimensions are assumed to be parameters
 * Destroys P
 * factor is NbEq x (nparam+2) matrix, containing stride constraints
 * on the parameters; column nparam is the constant;
 * column nparam+1 is the stride
 */
Polyhedron *remove_equalities_p(Polyhedron *P, unsigned nvar, Matrix **factor)
{
    Value g;
    Vector *v;
    Polyhedron *p = P, *q;
    unsigned dim = p->Dimension;
    Matrix *m1, *m2, *f;
    int i, j, skip;

    value_init(g);
    f = Matrix_Alloc(p->NbEq, dim-nvar+2);
    j = 0;
    *factor = f;
    skip = 0;
    while (nvar > 0 && p->NbEq - skip > 0) {
	assert(dim > 0);

	while (value_zero_p(p->Constraint[skip][0]) &&
	       First_Non_Zero(p->Constraint[skip]+1, nvar) == -1)
	    ++skip;
	if (p->NbEq == skip)
	    break;

	Vector_Gcd(p->Constraint[skip]+1, dim+1, &g);
	Vector_AntiScale(p->Constraint[skip]+1, p->Constraint[skip]+1, g, dim+1);
	Vector_Gcd(p->Constraint[skip]+1, nvar, &g);
	Vector_Copy(p->Constraint[skip]+1+nvar, f->p[j], dim-nvar+1);
	value_assign(f->p[j][dim-nvar+1], g);
	v = Vector_Alloc(dim);
	Vector_AntiScale(p->Constraint[skip]+1, v->p, g, nvar);
	Vector_Set(v->p+nvar, 0, dim-nvar);
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
	--nvar;
	++j;
    }
    value_clear(g);
    return p;
}

struct single {
    int	nr;
    int pos[2];
};

static void free_singles(int **singles, int dim)
{
    int i;
    for (i = 0; i < dim; ++i)
	free(singles[i]);
    free(singles);
}

static int **find_singles(Polyhedron *P, int dim, int max, int *nsingle)
{
    int i, j, prev;
    int **singles = (int **) malloc(dim * sizeof(int *));
    assert(singles);

    for (i = 0; i < dim; ++i) {
	singles[i] = (int *) malloc((max + 1) *sizeof(int));
	singles[i][0] = 0;
    }

    for (i = 0; i < P->NbConstraints; ++i) {
	for (j = 0, prev = -1; j < dim; ++j) {
	    if (value_notzero_p(P->Constraint[i][j+1])) {
		if (prev == -1)
		    prev = j;
		else {
		    if (prev != -2)
			singles[prev][0] = -1;
		    singles[j][0] = -1;
		    prev = -2;
		}
	    }
	}
	if (prev >= 0 && singles[prev][0] >= 0)
	    singles[prev][++singles[prev][0]] = i;
    }
    *nsingle = 0;
    for (j = 0; j < dim; ++j)
	if (singles[j][0] > 0)
	    ++*nsingle;
    if (!*nsingle) {
	free_singles(singles, dim);
	singles = 0;
    }
    return singles;
}

/*
 * The number of points in P is equal to factor time
 * the number of points in the polyhedron returned.
 * The return value is zero if no reduction can be found.
 */
Polyhedron* Polyhedron_Reduce(Polyhedron *P, Value* factor)
{
    int i, j, nsingle, k, p;
    unsigned dim = P->Dimension;
    int **singles;

    value_set_si(*factor, 1);

    assert (P->NbEq == 0);

    singles = find_singles(P, dim, 2, &nsingle);

    if (nsingle == 0)
	return 0;

    {
	Value tmp, pos, neg;
	Matrix *m = Matrix_Alloc((dim-nsingle)+1, dim+1);

	value_init(tmp);
	value_init(pos);
	value_init(neg);

	for (i = 0, j = 0; i < dim; ++i) {
	    if (singles[i][0] != 2)
		value_set_si(m->p[j++][i], 1);
	    else {
		for (k = 0; k <= 1; ++k) {
		    p = singles[i][1+k];
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
	free_singles(singles, dim);

	value_clear(tmp);
	value_clear(pos);
	value_clear(neg);
    }

    return P;
}

/*
 * Replaces constraint a x >= c by x >= ceil(c/a)
 * where "a" is a common factor in the coefficients
 * old is the constraint; v points to an initialized
 * value that this procedure can use.
 * Return non-zero if something changed.
 * Result is placed in new.
 */
int ConstraintSimplify(Value *old, Value *new, int len, Value* v)
{
    Vector_Gcd(old+1, len - 2, v);

    if (value_one_p(*v))
	return 0;

    Vector_AntiScale(old+1, new+1, *v, len-2);
    mpz_fdiv_q(new[len-1], old[len-1], *v);
    return 1;
}

/*
 * Project on final dim dimensions
 */
static Polyhedron* Polyhedron_Project(Polyhedron *P, int dim)
{
    int i;
    int remove = P->Dimension - dim;

    if (P->Dimension == dim)
	return Polyhedron_Copy(P);

    Matrix *T = Matrix_Alloc(dim+1, P->Dimension+1);
    for (i = 0; i < dim+1; ++i)
	value_set_si(T->p[i][i+remove], 1);
    Polyhedron *I = Polyhedron_Image(P, T, P->NbConstraints);
    Matrix_Free(T);
    return I;
}

struct section { Polyhedron * D; evalue E; };

static Polyhedron* ParamPolyhedron_Reduce_mod(Polyhedron *P, unsigned nvar, 
					      evalue* factor)
{
    int nsingle;
    int **singles;
    unsigned dim = P->Dimension;

    singles = find_singles(P, nvar, P->NbConstraints, &nsingle);

    if (nsingle == 0)
	return 0;

    {
	Polyhedron *C, *T;
	int i, j, p, n;
	Matrix *m = Matrix_Alloc((dim-nsingle)+1, dim+1);
	Value tmp, g;
	evalue mone;
	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
	C = Polyhedron_Project(P, dim-nvar);

	value_init(tmp);
	value_init(g);

	for (i = 0, j = 0; i < dim; ++i) {
	    if (i >= nvar || singles[i][0] < 2)
		value_set_si(m->p[j++][i], 1);
	    else {
		struct section *s;
		Matrix *M, *M2;
		int nd = 0;
		int k, l, k2, l2, q;
		evalue *L, *U;
		evalue F;
		/* put those with positive coefficients first; number: p */
		for (p = 0, n = singles[i][0]-1; p <= n; ) {
		    while (value_pos_p(P->Constraint[singles[i][1+p]][i+1]))
			++p;
		    while (value_neg_p(P->Constraint[singles[i][1+n]][i+1]))
			--n;
		    if (p < n) {
			int t = singles[i][1+p];
			singles[i][1+p] = singles[i][1+n];
			singles[i][1+n] = t;
			++p;
			--n;
		    }
		}
		n = singles[i][0]-p;
		assert (p >= 1 && n >= 1);
		s = (struct section *) malloc(p * n * sizeof(struct section));
		M = Matrix_Alloc((p-1) + (n-1), dim-nvar+2);
		for (k = 0; k < p; ++k) {
		    for (k2 = 0; k2 < p; ++k2) {
			if (k2 == k)
			    continue;
			q = k2 - (k2 > k);
			value_oppose(tmp, P->Constraint[singles[i][1+k2]][i+1]);
			value_set_si(M->p[q][0], 1);
			Vector_Combine(P->Constraint[singles[i][1+k]]+1+nvar,
				       P->Constraint[singles[i][1+k2]]+1+nvar,
				       M->p[q]+1,
				       tmp,
				       P->Constraint[singles[i][1+k]][i+1],
				       dim-nvar+1);
			if (k2 > k)
			    value_decrement(M->p[q][dim-nvar+1],
					    M->p[q][dim-nvar+1]);
			ConstraintSimplify(M->p[q], M->p[q], 
					   dim-nvar+2, &g);
		    }
		    for (l = p; l < p+n; ++l) {
			value_oppose(tmp, P->Constraint[singles[i][1+l]][i+1]);
			for (l2 = p; l2 < p+n; ++l2) {
			    if (l2 == l)
				continue;
			    q = l2-1 - (l2 > l);
			    value_set_si(M->p[q][0], 1);
			    Vector_Combine(P->Constraint[singles[i][1+l2]]+1+nvar,
					   P->Constraint[singles[i][1+l]]+1+nvar,
					   M->p[q]+1,
					   tmp,
					   P->Constraint[singles[i][1+l2]][i+1],
					   dim-nvar+1);
			    if (l2 > l)
				value_decrement(M->p[q][dim-nvar+1],
						M->p[q][dim-nvar+1]);
			    ConstraintSimplify(M->p[q], M->p[q], 
					       dim-nvar+2, &g);
			}
			M2 = Matrix_Copy(M);
			s[nd].D = Constraints2Polyhedron(M2, P->NbRays);
			Matrix_Free(M2);
			if (emptyQ(s[nd].D)) {
			    Polyhedron_Free(s[nd].D);
			    continue;
			}
			T = DomainIntersection(s[nd].D, C, 0);
			L = bv_ceil3(P->Constraint[singles[i][1+k]]+1+nvar, 
				  dim-nvar+1,
				  P->Constraint[singles[i][1+k]][i+1], T);
			U = bv_ceil3(P->Constraint[singles[i][1+l]]+1+nvar, 
				  dim-nvar+1,
				  P->Constraint[singles[i][1+l]][i+1], T);
			Domain_Free(T);
			eadd(L, U);
			eadd(&mone, U);
			emul(&mone, U);
			s[nd].E = *U;
			free_evalue_refs(L); 
			free(L);
			free(U);
			++nd;
		    }
		}

		Matrix_Free(M);

		value_init(F.d);
		value_set_si(F.d, 0);
		F.x.p = new_enode(partition, 2*nd, dim-nvar);
		for (k = 0; k < nd; ++k) {
		    EVALUE_SET_DOMAIN(F.x.p->arr[2*k], s[k].D);
		    value_clear(F.x.p->arr[2*k+1].d);
		    F.x.p->arr[2*k+1] = s[k].E;
		}
		free(s);

		emul(&F, factor);
		free_evalue_refs(&F); 
	    }
	}
	value_set_si(m->p[dim-nsingle][dim], 1);
	P = Polyhedron_Image(P, m, P->NbConstraints);
	Matrix_Free(m);
	free_singles(singles, nvar);

	value_clear(g);
	value_clear(tmp);

	free_evalue_refs(&mone); 
	Polyhedron_Free(C);
    }

    reduce_evalue(factor);

    return P;
}

#ifdef USE_MODULO
Polyhedron* ParamPolyhedron_Reduce(Polyhedron *P, unsigned nvar, 
				   evalue* factor)
{
    return ParamPolyhedron_Reduce_mod(P, nvar, factor);
}
#else
Polyhedron* ParamPolyhedron_Reduce(Polyhedron *P, unsigned nvar, 
				   evalue* factor)
{
    Polyhedron *R;
    evalue tmp;
    value_init(tmp.d);
    evalue_set_si(&tmp, 1, 1);
    R = ParamPolyhedron_Reduce_mod(P, nvar, &tmp);
    evalue_mod2table(&tmp, P->Dimension - nvar);
    reduce_evalue(&tmp);
    emul(&tmp, factor);
    free_evalue_refs(&tmp); 
    return R;
}
#endif

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

void Enumeration_Print(FILE *Dst, Enumeration *en, char **params)
{
    for (; en; en = en->next) {
	Print_Domain(Dst, en->ValidityDomain, params);
	print_evalue(Dst, &en->EP, params);
    }
}

void Enumeration_mod2table(Enumeration *en, unsigned nparam)
{
    for (; en; en = en->next) {
	evalue_mod2table(&en->EP, nparam);
	reduce_evalue(&en->EP);
    }
}

size_t Enumeration_size(Enumeration *en)
{
    size_t s = 0;

    for (; en; en = en->next) {
	s += domain_size(en->ValidityDomain);
	s += evalue_size(&en->EP);
    }
    return s;
}

void Free_ParamNames(char **params, int m)
{
    while (--m >= 0)
	free(params[m]);
    free(params);
}

int DomainIncludes(Polyhedron *Pol1, Polyhedron *Pol2)
{
    Polyhedron *P2;
    for ( ; Pol1; Pol1 = Pol1->next) {
	for (P2 = Pol2; P2; P2 = P2->next)
	    if (!PolyhedronIncludes(Pol1, P2))
		break;
	if (!P2)
	    return 1;
    }
    return 0;
}

static Polyhedron *p_simplify_constraints(Polyhedron *P, Vector *row,
					  Value *g, unsigned MaxRays)
{
    Polyhedron *T, *R = P;
    int len = P->Dimension+2;
    int r;

    /* Also look at equalities.
     * If an equality can be "simplified" then there
     * are no integer solutions anyway and the following loop
     * will add a conflicting constraint
     */
    for (r = 0; r < R->NbConstraints; ++r) {
	if (ConstraintSimplify(R->Constraint[r], row->p, len, g)) {
	    T = R;
	    R = AddConstraints(row->p, 1, R, MaxRays);
	    if (T != P)
		Polyhedron_Free(T);
	    r = -1;
	}
    }
    if (R != P)
	Polyhedron_Free(P);
    return R;
}

/*
 * Replaces constraint a x >= c by x >= ceil(c/a)
 * where "a" is a common factor in the coefficients
 * Destroys P and returns a newly allocated Polyhedron
 * or just returns P in case no changes were made
 */
Polyhedron *DomainConstraintSimplify(Polyhedron *P, unsigned MaxRays)
{
    Polyhedron **prev;
    int len = P->Dimension+2;
    Vector *row = Vector_Alloc(len);
    value_set_si(row->p[0], 1);
    Value g;
    value_init(g);
    Polyhedron *R = P, *N;

    for (prev = &R; P; P = N) {
	Polyhedron *T;
	N = P->next;
	T = p_simplify_constraints(P, row, &g, MaxRays);

	if (emptyQ(T) && prev != &R) {
	    Polyhedron_Free(T);
	    *prev = NULL;
	    continue;
	}

	if (T != P)
	    T->next = N;
	*prev = T;
	prev = &T->next;
    }

    if (R->next && emptyQ(R)) {
	N = R->next;
	Polyhedron_Free(R);
	R = N;
    }

    value_clear(g);
    Vector_Free(row);
    return R;
}

int line_minmax(Polyhedron *I, Value *min, Value *max)
{
    int i;

    if (I->NbEq >= 1) {
	value_oppose(I->Constraint[0][2], I->Constraint[0][2]);
	/* There should never be a remainder here */
	if (value_pos_p(I->Constraint[0][1]))
	    mpz_fdiv_q(*min, I->Constraint[0][2], I->Constraint[0][1]);
	else
	    mpz_fdiv_q(*min, I->Constraint[0][2], I->Constraint[0][1]);
	value_assign(*max, *min);
    } else for (i = 0; i < I->NbConstraints; ++i) {
	if (value_zero_p(I->Constraint[i][1])) {
	    Polyhedron_Free(I);
	    return 0;
	}

	value_oppose(I->Constraint[i][2], I->Constraint[i][2]);
	if (value_pos_p(I->Constraint[i][1]))
	    mpz_cdiv_q(*min, I->Constraint[i][2], I->Constraint[i][1]);
	else
	    mpz_fdiv_q(*max, I->Constraint[i][2], I->Constraint[i][1]);
    }
    Polyhedron_Free(I);
    return 1;
}

/** 

PROCEDURES TO COMPUTE ENUMERATION. recursive procedure, recurse for
each imbriquation

@param pos index position of current loop index (1..hdim-1)
@param P loop domain
@param context context values for fixed indices 
@param exist	number of existential variables
@return the number of integer points in this
polyhedron 

*/
void count_points_e (int pos, Polyhedron *P, int exist, int nparam,
		     Value *context, Value *res)
{
    Value LB, UB, k, c;

    if (emptyQ(P)) {
	value_set_si(*res, 0);
	return;
    }

    value_init(LB); value_init(UB); value_init(k);
    value_set_si(LB,0);
    value_set_si(UB,0);

    if (lower_upper_bounds(pos,P,context,&LB,&UB) !=0) {
        /* Problem if UB or LB is INFINITY */
        value_clear(LB); value_clear(UB); value_clear(k);
	if (pos > P->Dimension - nparam - exist)
	    value_set_si(*res, 1);
	else 
	    value_set_si(*res, -1);
	return;
    }
  
#ifdef EDEBUG1
    if (!P->next) {
        int i;
        for (value_assign(k,LB); value_le(k,UB); value_increment(k,k)) {
            fprintf(stderr, "(");
            for (i=1; i<pos; i++) {
                value_print(stderr,P_VALUE_FMT,context[i]);
                fprintf(stderr,",");
            }
            value_print(stderr,P_VALUE_FMT,k);
            fprintf(stderr,")\n");
        }
    }
#endif
  
    value_set_si(context[pos],0);
    if (value_lt(UB,LB)) {
        value_clear(LB); value_clear(UB); value_clear(k);
	value_set_si(*res, 0);
	return;  
    }  
    if (!P->next) {
	if (exist)
	    value_set_si(*res, 1);
	else {
	    value_substract(k,UB,LB);
	    value_add_int(k,k,1);
	    value_assign(*res, k);
	}
        value_clear(LB); value_clear(UB); value_clear(k);
        return;
    } 

    /*-----------------------------------------------------------------*/
    /* Optimization idea                                               */
    /*   If inner loops are not a function of k (the current index)    */
    /*   i.e. if P->Constraint[i][pos]==0 for all P following this and */
    /*        for all i,                                               */
    /*   Then CNT = (UB-LB+1)*count_points(pos+1, P->next, context)    */
    /*   (skip the for loop)                                           */
    /*-----------------------------------------------------------------*/
  
    value_init(c);
    value_set_si(*res, 0);
    for (value_assign(k,LB);value_le(k,UB);value_increment(k,k)) {
        /* Insert k in context */
        value_assign(context[pos],k);
	count_points_e(pos+1, P->next, exist, nparam, context, &c);
	if(value_notmone_p(c))
	    value_addto(*res, *res, c);
	else {
	    value_set_si(*res, -1);
	    break;
        }
	if (pos > P->Dimension - nparam - exist &&
		value_pos_p(*res))
	    break;
    }
    value_clear(c);
  
#ifdef EDEBUG11
    fprintf(stderr,"%d\n",CNT);
#endif
  
    /* Reset context */
    value_set_si(context[pos],0);
    value_clear(LB); value_clear(UB); value_clear(k);
    return;
} /* count_points_e */
