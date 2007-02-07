#include <assert.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <barvinok/barvinok.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include "config.h"
#include "conversion.h"
#include "lattice_point.h"

using std::cerr;
using std::endl;

#define ALLOC(type) (type*)malloc(sizeof(type))

/* returns an evalue that corresponds to
 *
 * c/(*den) x_param
 */
static evalue *term(int param, ZZ& c, Value *den = NULL)
{
    evalue *EP = new evalue();
    value_init(EP->d);
    value_set_si(EP->d,0);
    EP->x.p = new_enode(polynomial, 2, param + 1);
    evalue_set_si(&EP->x.p->arr[0], 0, 1);
    value_init(EP->x.p->arr[1].x.n);
    if (den == NULL)
	value_set_si(EP->x.p->arr[1].d, 1);
    else
	value_assign(EP->x.p->arr[1].d, *den);
    zz2value(c, EP->x.p->arr[1].x.n);
    return EP;
}

/* returns an evalue that corresponds to
 *
 *   sum_i p[i] * x_i
 */
evalue *multi_monom(vec_ZZ& p)
{
    evalue *X = new evalue();
    value_init(X->d);
    value_init(X->x.n);
    unsigned nparam = p.length()-1;
    zz2value(p[nparam], X->x.n);
    value_set_si(X->d, 1);
    for (int i = 0; i < nparam; ++i) {
	if (p[i] == 0)
	    continue;
	evalue *T = term(i, p[i]);
	eadd(T, X); 
	free_evalue_refs(T); 
	delete T;
    }
    return X;
}

/*
 * Check whether mapping polyhedron P on the affine combination
 * num yields a range that has a fixed quotient on integer
 * division by d
 * If zero is true, then we are only interested in the quotient
 * for the cases where the remainder is zero.
 * Returns NULL if false and a newly allocated value if true.
 */
static Value *fixed_quotient(Polyhedron *P, vec_ZZ& num, Value d, bool zero)
{
    Value* ret = NULL;
    int len = num.length();
    Matrix *T = Matrix_Alloc(2, len);
    zz2values(num, T->p[0]);
    value_set_si(T->p[1][len-1], 1);
    Polyhedron *I = Polyhedron_Image(P, T, P->NbConstraints);
    Matrix_Free(T);

    int i;
    for (i = 0; i < I->NbRays; ++i)
	if (value_zero_p(I->Ray[i][2])) {
	    Polyhedron_Free(I);
	    return NULL;
	}

    Value min, max;
    value_init(min);
    value_init(max);
    int bounded = line_minmax(I, &min, &max);
    assert(bounded);

    if (zero)
	mpz_cdiv_q(min, min, d);
    else
	mpz_fdiv_q(min, min, d);
    mpz_fdiv_q(max, max, d);

    if (value_eq(min, max)) {
	ret = ALLOC(Value);
	value_init(*ret);
	value_assign(*ret, min);
    } 
    value_clear(min);
    value_clear(max);
    return ret;
}

/*
 * Normalize linear expression coef modulo m
 * Removes common factor and reduces coefficients
 * Returns index of first non-zero coefficient or len
 */
int normal_mod(Value *coef, int len, Value *m)
{
    Value gcd;
    value_init(gcd);

    Vector_Gcd(coef, len, &gcd);
    Gcd(gcd, *m, &gcd);
    Vector_AntiScale(coef, coef, gcd, len);

    value_division(*m, *m, gcd);
    value_clear(gcd);

    if (value_one_p(*m))
	return len;

    int j;
    for (j = 0; j < len; ++j)
	mpz_fdiv_r(coef[j], coef[j], *m);
    for (j = 0; j < len; ++j)
	if (value_notzero_p(coef[j]))
	    break;

    return j;
}

static bool mod_needed(Polyhedron *PD, vec_ZZ& num, Value d, evalue *E)
{
    Value *q = fixed_quotient(PD, num, d, false);

    if (!q)
	return true;

    value_oppose(*q, *q);
    evalue EV;
    value_init(EV.d);
    value_set_si(EV.d, 1);
    value_init(EV.x.n);
    value_multiply(EV.x.n, *q, d);
    eadd(&EV, E);
    free_evalue_refs(&EV); 
    value_clear(*q);
    free(q);
    return false;
}

/* modifies f argument ! */
static void ceil_mod(Value *coef, int len, Value d, ZZ& f, evalue *EP, Polyhedron *PD)
{
    Value m;
    value_init(m);
    value_set_si(m, -1);

    Vector_Scale(coef, coef, m, len);

    value_assign(m, d);
    int j = normal_mod(coef, len, &m);

    if (j == len) {
	value_clear(m);
	return;
    }

    vec_ZZ num;
    values2zz(coef, num, len);

    ZZ g;
    value2zz(m, g);

    evalue tmp;
    value_init(tmp.d);
    evalue_set_si(&tmp, 0, 1);

    int p = j;
    if (g % 2 == 0)
	while (j < len-1 && (num[j] == g/2 || num[j] == 0))
	    ++j;
    if ((j < len-1 && num[j] > g/2) || (j == len-1 && num[j] >= (g+1)/2)) {
	for (int k = j; k < len-1; ++k)
	    if (num[k] != 0)
		num[k] = g - num[k];
	num[len-1] = g - 1 - num[len-1];
	value_assign(tmp.d, m);
	ZZ t = f*(g-1);
	zz2value(t, tmp.x.n);
	eadd(&tmp, EP);
	f = -f;
    }

    if (p >= len-1) {
	ZZ t = num[len-1] * f;
	zz2value(t, tmp.x.n);
	value_assign(tmp.d, m);
	eadd(&tmp, EP);
    } else {
	evalue *E = multi_monom(num);
	evalue EV;
	value_init(EV.d);

	if (PD && !mod_needed(PD, num, m, E)) {
	    value_init(EV.x.n);
	    zz2value(f, EV.x.n);
	    value_assign(EV.d, m);
	    emul(&EV, E);
	    eadd(E, EP);
	} else {
	    value_init(EV.x.n);
	    value_set_si(EV.x.n, 1);
	    value_assign(EV.d, m);
	    emul(&EV, E);
	    value_clear(EV.x.n);
	    value_set_si(EV.d, 0);
	    EV.x.p = new_enode(fractional, 3, -1);
	    evalue_copy(&EV.x.p->arr[0], E);
	    evalue_set_si(&EV.x.p->arr[1], 0, 1);
	    value_init(EV.x.p->arr[2].x.n);
	    zz2value(f, EV.x.p->arr[2].x.n);
	    value_set_si(EV.x.p->arr[2].d, 1);

	    eadd(&EV, EP);
	}

	free_evalue_refs(&EV); 
	free_evalue_refs(E); 
	delete E;
    }

    free_evalue_refs(&tmp); 

out:
    value_clear(m);
}

static void ceil(Value *coef, int len, Value d, ZZ& f,
                 evalue *EP, Polyhedron *PD, barvinok_options *options)
{
    ceil_mod(coef, len, d, f, EP, PD);
    if (options->lookup_table)
	evalue_mod2table(EP, len-1);
}

evalue* bv_ceil3(Value *coef, int len, Value d, Polyhedron *P)
{
    Vector *val = Vector_Alloc(len);

    Value t;
    value_init(t);
    value_set_si(t, -1);
    Vector_Scale(coef, val->p, t, len);
    value_absolute(t, d);

    vec_ZZ num;
    values2zz(val->p, num, len);
    evalue *EP = multi_monom(num);

    evalue tmp;
    value_init(tmp.d);
    value_init(tmp.x.n);
    value_set_si(tmp.x.n, 1);
    value_assign(tmp.d, t);

    emul(&tmp, EP);

    ZZ one;
    one = 1;
    ceil_mod(val->p, len, t, one, EP, P);
    value_clear(t);

    /* copy EP to malloc'ed evalue */
    evalue *E = ALLOC(evalue);
    *E = *EP;
    delete EP;

    free_evalue_refs(&tmp); 
    Vector_Free(val);

    return E;
}

void lattice_point(Value* values, const mat_ZZ& rays, vec_ZZ& vertex, int *closed)
{
    unsigned dim = rays.NumRows();
    if (value_one_p(values[dim]) && !closed)
	values2zz(values, vertex, dim);
    else {
	Matrix* Rays = rays2matrix(rays);
	Matrix *inv = Matrix_Alloc(Rays->NbRows, Rays->NbColumns);
	int ok = Matrix_Inverse(Rays, inv);
	assert(ok);
	Matrix_Free(Rays);
	Rays = rays2matrix(rays);
	Vector *lambda = Vector_Alloc(dim+1);
	Vector_Matrix_Product(values, inv, lambda->p);
	Matrix_Free(inv);
	for (int j = 0; j < dim; ++j)
	    if (!closed || closed[j])
		mpz_cdiv_q(lambda->p[j], lambda->p[j], lambda->p[dim]);
	    else {
		value_addto(lambda->p[j], lambda->p[j], lambda->p[dim]);
		mpz_fdiv_q(lambda->p[j], lambda->p[j], lambda->p[dim]);
	    }
	value_set_si(lambda->p[dim], 1);
	Vector *A = Vector_Alloc(dim+1);
	Vector_Matrix_Product(lambda->p, Rays, A->p);
	Vector_Free(lambda);
	Matrix_Free(Rays);
	values2zz(A->p, vertex, dim);
	Vector_Free(A);
    }
}

/* Compute the lattice points in the vertex cone at "values" with rays "rays".
 * The lattice points are returned in "vertex".
 *
 * Rays has the generators as rows and so does W.
 * We first compute { m-v, u_i^* } with m = k W, where k runs through
 * the cosets.
 * We compute
 * [k 1] [ d1*W  0 ] [  U'  0 ] = [k 1] T2
 *       [ -v   d1 ] [  0  d2 ]
 * where d1 and d2 are the denominators of v and U^{-1}=U'/d2.
 * Then  lambda = { k } (componentwise)
 * We compute x - floor(x) = {x} = { a/b } as fdiv_r(a,b)/b
 * For open rays/facets, we need values in (0,1] rather than [0,1),
 * so we compute {{x}} = x - ceil(x-1) = a/b - ceil((a-b)/b)
 *		       = (a - b cdiv_q(a-b,b) - b + b)/b
 *		       = (cdiv_r(a,b)+b)/b
 * Finally, we compute v + lambda * U
 * The denominator of lambda can be d1*d2, that of lambda2 = lambda*U
 * can be at most d1, since it is integer if v = 0.
 * The denominator of v + lambda2 is 1.
 */
void lattice_point(Value* values, const mat_ZZ& rays, mat_ZZ& vertex,
		   unsigned long det, int *closed)
{
    unsigned dim = rays.NumRows();
    vertex.SetDims(det, dim);
    if (det == 1) {
	lattice_point(values, rays, vertex[0], closed);
	return;
    }
    Matrix* Rays = rays2matrix2(rays);
    Matrix *U, *W, *D;
    Smith(Rays, &U, &W, &D);
    Matrix_Free(Rays);
    Matrix_Free(U);

    Matrix *T = Matrix_Alloc(W->NbRows+1, W->NbColumns+1);
    for (int i = 0; i < W->NbRows; ++i)
	Vector_Scale(W->p[i], T->p[i], values[dim], W->NbColumns);
    Matrix_Free(W);
    Value tmp;
    value_init(tmp);
    value_set_si(tmp, -1);
    Vector_Scale(values, T->p[dim], tmp, dim);
    value_clear(tmp);
    value_assign(T->p[dim][dim], values[dim]);

    Rays = rays2matrix(rays);
    Matrix *inv = Matrix_Alloc(Rays->NbRows, Rays->NbColumns);
    int ok = Matrix_Inverse(Rays, inv);
    assert(ok);
    Matrix_Free(Rays);

    Matrix *T2 = Matrix_Alloc(dim+1, dim+1);
    Matrix_Product(T, inv, T2);
    Matrix_Free(T);

    Rays = rays2matrix(rays);

    Vector *k = Vector_Alloc(dim+1);
    value_set_si(k->p[dim], 1);
    Vector *lambda = Vector_Alloc(dim+1);
    Vector *lambda2 = Vector_Alloc(dim+1);
    for (unsigned long i = 0; i < det; ++i) {
	unsigned long val = i;
	for (int j = 0; j < dim; ++j) {
	    value_set_si(k->p[j], val % mpz_get_ui(D->p[j][j]));
	    val /= mpz_get_ui(D->p[j][j]);
	}
	Vector_Matrix_Product(k->p, T2, lambda->p);
	for (int j = 0; j < dim; ++j)
	    if (!closed || closed[j])
		mpz_fdiv_r(lambda->p[j], lambda->p[j], lambda->p[dim]);
	    else {
		mpz_cdiv_r(lambda->p[j], lambda->p[j], lambda->p[dim]);
		value_addto(lambda->p[j], lambda->p[j], lambda->p[dim]);
	    }
	Vector_Matrix_Product(lambda->p, Rays, lambda2->p);
	for (int j = 0; j < dim; ++j)
	    assert(mpz_divisible_p(lambda2->p[j], inv->p[dim][dim]));
	Vector_AntiScale(lambda2->p, lambda2->p, inv->p[dim][dim], dim+1);
	Vector_Add(lambda2->p, values, lambda2->p, dim);
	for (int j = 0; j < dim; ++j)
	    assert(mpz_divisible_p(lambda2->p[j], values[dim]));
	Vector_AntiScale(lambda2->p, lambda2->p, values[dim], dim+1);
	values2zz(lambda2->p, vertex[i], dim);
    }
    Vector_Free(k);
    Vector_Free(lambda);
    Vector_Free(lambda2);
    Matrix_Free(D);
    Matrix_Free(Rays);
    Matrix_Free(inv);

    Matrix_Free(T2);
}

static void vertex_period(
		    const mat_ZZ& rays, vec_ZZ& lambda, Matrix *T, 
		    Value lcm, int p, Vector *val, 
		    evalue *E, evalue* ev,
		    ZZ& offset)
{
    unsigned nparam = T->NbRows - 1;
    unsigned dim = rays.NumRows();
    Value tmp;
    ZZ nump;

    if (p == nparam) {
	vec_ZZ vertex;
	ZZ num, l;
	Vector * values = Vector_Alloc(dim + 1);
	Vector_Matrix_Product(val->p, T, values->p);
	value_assign(values->p[dim], lcm);
	lattice_point(values->p, rays, vertex, NULL);
	num = vertex * lambda;
	value2zz(lcm, l);
	num *= l;
	num += offset;
	value_init(ev->x.n);
	zz2value(num, ev->x.n);
	value_assign(ev->d, lcm);
	Vector_Free(values);
	return;
    }

    value_init(tmp);
    vec_ZZ vertex;
    values2zz(T->p[p], vertex, dim);
    nump = vertex * lambda;
    if (First_Non_Zero(val->p, p) == -1) {
	value_assign(tmp, lcm);
	evalue *ET = term(p, nump, &tmp);
	eadd(ET, E);   
	free_evalue_refs(ET); 
	delete ET;
    }

    value_assign(tmp, lcm);
    if (First_Non_Zero(T->p[p], dim) != -1)
	Vector_Gcd(T->p[p], dim, &tmp);
    Gcd(tmp, lcm, &tmp);
    if (value_lt(tmp, lcm)) {
	ZZ count;

	value_division(tmp, lcm, tmp);
	value_set_si(ev->d, 0);
	ev->x.p = new_enode(periodic, VALUE_TO_INT(tmp), p+1);
	value2zz(tmp, count);
	do {
	    value_decrement(tmp, tmp);
	    --count;
	    ZZ new_offset = offset - count * nump;
	    value_assign(val->p[p], tmp);
	    vertex_period(rays, lambda, T, lcm, p+1, val, E, 
			  &ev->x.p->arr[VALUE_TO_INT(tmp)], new_offset);
	} while (value_pos_p(tmp));
    } else
	vertex_period(rays, lambda, T, lcm, p+1, val, E, ev, offset);
    value_clear(tmp);
}

/* Returns the power of (t+1) in the term of a rational generating function,
 * i.e., the scalar product of the actual lattice point and lambda.
 * The lattice point is the unique lattice point in the fundamental parallelepiped
 * of the unimodual cone i shifted to the parametric vertex W/lcm.
 *
 * The rows of W refer to the coordinates of the vertex
 * The first nparam columns are the coefficients of the parameters
 * and the final column is the constant term.
 * lcm is the common denominator of all coefficients.
 *
 * PD is the parameter domain, which, if != NULL, may be used to simply the
 * resulting expression.
 */
static evalue* lattice_point_fractional(const mat_ZZ& rays, vec_ZZ& lambda,
					Matrix *W, Value lcm, Polyhedron *PD)
{
    unsigned nparam = W->NbColumns - 1;

    Matrix* Rays = rays2matrix2(rays);
    Matrix *T = Transpose(Rays);
    Matrix *T2 = Matrix_Copy(T);
    Matrix *inv = Matrix_Alloc(T2->NbRows, T2->NbColumns);
    int ok = Matrix_Inverse(T2, inv);
    assert(ok);
    Matrix_Free(Rays);
    Matrix_Free(T2);
    mat_ZZ vertex;
    matrix2zz(W, vertex, W->NbRows, W->NbColumns);

    vec_ZZ num;
    num = lambda * vertex;

    evalue *EP = multi_monom(num);

    evalue tmp;
    value_init(tmp.d);
    value_init(tmp.x.n);
    value_set_si(tmp.x.n, 1);
    value_assign(tmp.d, lcm);

    emul(&tmp, EP);

    Matrix *L = Matrix_Alloc(inv->NbRows, W->NbColumns);
    Matrix_Product(inv, W, L);

    mat_ZZ RT;
    matrix2zz(T, RT, T->NbRows, T->NbColumns);
    Matrix_Free(T);

    vec_ZZ p = lambda * RT;

    for (int i = 0; i < L->NbRows; ++i) {
	ceil_mod(L->p[i], nparam+1, lcm, p[i], EP, PD);
    }

    Matrix_Free(L);

    Matrix_Free(inv);
    free_evalue_refs(&tmp); 
    return EP;
}

static evalue* lattice_point_table(const mat_ZZ& rays, vec_ZZ& lambda, Matrix *W,
				   Value lcm, Polyhedron *PD)
{
    Matrix *T = Transpose(W);
    unsigned nparam = T->NbRows - 1;

    evalue *EP = new evalue();
    value_init(EP->d);
    evalue_set_si(EP, 0, 1);

    evalue ev;
    Vector *val = Vector_Alloc(nparam+1);
    value_set_si(val->p[nparam], 1);
    ZZ offset(INIT_VAL, 0);
    value_init(ev.d);
    vertex_period(rays, lambda, T, lcm, 0, val, EP, &ev, offset);
    Vector_Free(val);
    eadd(&ev, EP);
    free_evalue_refs(&ev);   

    Matrix_Free(T);

    reduce_evalue(EP);

    return EP;
}

evalue* lattice_point(const mat_ZZ& rays, vec_ZZ& lambda, Matrix *W,
		      Value lcm, Polyhedron *PD, barvinok_options *options)
{
    if (options->lookup_table)
	return lattice_point_table(rays, lambda, W, lcm, PD);
    else
	return lattice_point_fractional(rays, lambda, W, lcm, PD);
}

/* returns the unique lattice point in the fundamental parallelepiped
 * of the unimodual cone C shifted to the parametric vertex V.
 *
 * The return values num and E_vertex are such that
 * coordinate i of this lattice point is equal to
 *
 *	    num[i] + E_vertex[i]
 */
void lattice_point(Param_Vertices *V, const mat_ZZ& rays, vec_ZZ& num, 
		   evalue **E_vertex, barvinok_options *options)
{
    unsigned nparam = V->Vertex->NbColumns - 2;
    unsigned dim = rays.NumCols();
    vec_ZZ vertex;
    vertex.SetLength(nparam+1);

    Value lcm, tmp;
    value_init(lcm);
    value_init(tmp);
    value_set_si(lcm, 1);

    for (int j = 0; j < V->Vertex->NbRows; ++j) {
	value_lcm(lcm, V->Vertex->p[j][nparam+1], &lcm);
    }

    if (value_notone_p(lcm)) {
	Matrix * mv = Matrix_Alloc(dim, nparam+1);
	for (int j = 0 ; j < dim; ++j) {
	    value_division(tmp, lcm, V->Vertex->p[j][nparam+1]);
	    Vector_Scale(V->Vertex->p[j], mv->p[j], tmp, nparam+1);
	}

	Matrix* Rays = rays2matrix2(rays);
	Matrix *T = Transpose(Rays);
	Matrix *T2 = Matrix_Copy(T);
	Matrix *inv = Matrix_Alloc(T2->NbRows, T2->NbColumns);
	int ok = Matrix_Inverse(T2, inv);
	assert(ok);
	Matrix_Free(Rays);
	Matrix_Free(T2);
	Matrix *L = Matrix_Alloc(inv->NbRows, mv->NbColumns);
	Matrix_Product(inv, mv, L);
	Matrix_Free(inv);

	evalue f;
	value_init(f.d);
	value_init(f.x.n);

	ZZ one;

	evalue *remainders[dim];
	for (int i = 0; i < dim; ++i) {
	    remainders[i] = evalue_zero();
	    one = 1;
	    ceil(L->p[i], nparam+1, lcm, one, remainders[i], 0, options);
	}
	Matrix_Free(L);


	for (int i = 0; i < V->Vertex->NbRows; ++i) {
	    values2zz(mv->p[i], vertex, nparam+1);
	    E_vertex[i] = multi_monom(vertex);
	    num[i] = 0;

	    value_set_si(f.x.n, 1);
	    value_assign(f.d, lcm);

	    emul(&f, E_vertex[i]);

	    for (int j = 0; j < dim; ++j) {
		if (value_zero_p(T->p[i][j]))
		    continue;
		evalue cp;
		value_init(cp.d);
		evalue_copy(&cp, remainders[j]);
		if (value_notone_p(T->p[i][j])) {
		    value_set_si(f.d, 1);
		    value_assign(f.x.n, T->p[i][j]);
		    emul(&f, &cp);
		}
		eadd(&cp, E_vertex[i]);
		free_evalue_refs(&cp);
	    }
	}
	for (int i = 0; i < dim; ++i) {
	    free_evalue_refs(remainders[i]); 
	    free(remainders[i]);
	}

	free_evalue_refs(&f); 

	Matrix_Free(T);
	Matrix_Free(mv);
	value_clear(lcm);
	value_clear(tmp);
	return;
    }
    value_clear(lcm);
    value_clear(tmp);

    for (int i = 0; i < V->Vertex->NbRows; ++i) {
	/* fixed value */
	if (First_Non_Zero(V->Vertex->p[i], nparam) == -1) {
	    E_vertex[i] = 0;
	    value2zz(V->Vertex->p[i][nparam], num[i]);
	} else {
	    values2zz(V->Vertex->p[i], vertex, nparam+1);
	    E_vertex[i] = multi_monom(vertex);
	    num[i] = 0;
	}
    }
}

