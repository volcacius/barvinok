#include <assert.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <barvinok/util.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>
#include "piputil.h"
}
#include "config.h"
#include <barvinok/barvinok.h>
#include <barvinok/genfun.h>
#include "conversion.h"
#include "decomposer.h"
#include "lattice_point.h"
#include "reduce_domain.h"
#include "genfun_constructor.h"
#include "sample.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::deque;
using std::string;
using std::ostringstream;

#define ALLOC(t,p) p = (t*)malloc(sizeof(*p))

static void rays(mat_ZZ& r, Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);
    r.SetDims(dim, dim);
    ZZ tmp;

    int i, c;
    for (i = 0, c = 0; i < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1])) {
	    for (int j = 0; j < dim; ++j) {
		value2zz(C->Ray[i][j+1], tmp);
		r[j][c] = tmp;
	    }
	    ++c;
	}
}

class dpoly_n {
public:
    Matrix *coeff;
    ~dpoly_n() {
	Matrix_Free(coeff);
    }
    dpoly_n(int d, ZZ& degree_0, ZZ& degree_1, int offset = 0) {
	Value d0, d1;
	value_init(d0);
	value_init(d1);
	zz2value(degree_0, d0);
	zz2value(degree_1, d1);
	coeff = Matrix_Alloc(d+1, d+1+1);
	value_set_si(coeff->p[0][0], 1);
	value_set_si(coeff->p[0][d+1], 1);
	for (int i = 1; i <= d; ++i) {
	    value_multiply(coeff->p[i][0], coeff->p[i-1][0], d0);
	    Vector_Combine(coeff->p[i-1], coeff->p[i-1]+1, coeff->p[i]+1,
			   d1, d0, i);
	    value_set_si(coeff->p[i][d+1], i);
	    value_multiply(coeff->p[i][d+1], coeff->p[i][d+1], coeff->p[i-1][d+1]);
	    value_decrement(d0, d0);
	}
	value_clear(d0);
	value_clear(d1);
    }
    void div(dpoly& d, Vector *count, ZZ& sign) {
	int len = coeff->NbRows;
	Matrix * c = Matrix_Alloc(coeff->NbRows, coeff->NbColumns);
	Value tmp;
	value_init(tmp);
	for (int i = 0; i < len; ++i) {
	    Vector_Copy(coeff->p[i], c->p[i], len+1);
	    for (int j = 1; j <= i; ++j) {
		zz2value(d.coeff[j], tmp);
		value_multiply(tmp, tmp, c->p[i][len]);
		value_oppose(tmp, tmp);
		Vector_Combine(c->p[i], c->p[i-j], c->p[i],
			       c->p[i-j][len], tmp, len);
		value_multiply(c->p[i][len], c->p[i][len], c->p[i-j][len]);
	    }
	    zz2value(d.coeff[0], tmp);
	    value_multiply(c->p[i][len], c->p[i][len], tmp);
	}
	if (sign == -1) {
	    value_set_si(tmp, -1);
	    Vector_Scale(c->p[len-1], count->p, tmp, len);
	    value_assign(count->p[len], c->p[len-1][len]);
	} else
	    Vector_Copy(c->p[len-1], count->p, len+1);
	Vector_Normalize(count->p, len+1);
	value_clear(tmp);
	Matrix_Free(c);
    }
};

const int MAX_TRY=10;
/*
 * Searches for a vector that is not orthogonal to any
 * of the rays in rays.
 */
static void nonorthog(mat_ZZ& rays, vec_ZZ& lambda)
{
    int dim = rays.NumCols();
    bool found = false;
    lambda.SetLength(dim);
    if (dim == 0)
	return;

    for (int i = 2; !found && i <= 50*dim; i+=4) {
	for (int j = 0; j < MAX_TRY; ++j) {
	    for (int k = 0; k < dim; ++k) {
		int r = random_int(i)+2;
		int v = (2*(r%2)-1) * (r >> 1);
		lambda[k] = v;
	    }
	    int k = 0;
	    for (; k < rays.NumRows(); ++k)
		if (lambda * rays[k] == 0)
		    break;
	    if (k == rays.NumRows()) {
		found = true;
		break;
	    }
	}
    }
    assert(found);
}

static void add_rays(mat_ZZ& rays, Polyhedron *i, int *r, int nvar = -1, 
		     bool all = false)
{
    unsigned dim = i->Dimension;
    if (nvar == -1)
	nvar = dim;
    for (int k = 0; k < i->NbRays; ++k) {
	if (!value_zero_p(i->Ray[k][dim+1]))
	    continue;
	if (!all && nvar != dim && First_Non_Zero(i->Ray[k]+1, nvar) == -1)
	    continue;
	values2zz(i->Ray[k]+1, rays[(*r)++], nvar);
    }
}

static void mask_r(Matrix *f, int nr, Vector *lcm, int p, Vector *val, evalue *ev)
{
    unsigned nparam = lcm->Size;

    if (p == nparam) {
	Vector * prod = Vector_Alloc(f->NbRows);
	Matrix_Vector_Product(f, val->p, prod->p);
	int isint = 1;
	for (int i = 0; i < nr; ++i) {
	    value_modulus(prod->p[i], prod->p[i], f->p[i][nparam+1]);
	    isint &= value_zero_p(prod->p[i]);
	}
	value_set_si(ev->d, 1);
	value_init(ev->x.n);
	value_set_si(ev->x.n, isint);
	Vector_Free(prod);
	return;
    }

    Value tmp;
    value_init(tmp);
    if (value_one_p(lcm->p[p]))
	mask_r(f, nr, lcm, p+1, val, ev);
    else { 
	value_assign(tmp, lcm->p[p]);
	value_set_si(ev->d, 0);
	ev->x.p = new_enode(periodic, VALUE_TO_INT(tmp), p+1);
	do {
	    value_decrement(tmp, tmp);
	    value_assign(val->p[p], tmp);
	    mask_r(f, nr, lcm, p+1, val, &ev->x.p->arr[VALUE_TO_INT(tmp)]);
	} while (value_pos_p(tmp));
    }
    value_clear(tmp);
}

#ifdef USE_MODULO
static void mask(Matrix *f, evalue *factor)
{
    int nr = f->NbRows, nc = f->NbColumns;
    int n;
    bool found = false;
    for (n = 0; n < nr && value_notzero_p(f->p[n][nc-1]); ++n)
	if (value_notone_p(f->p[n][nc-1]) &&
	    value_notmone_p(f->p[n][nc-1]))
		found = true;
    if (!found)
	return;

    evalue EP;
    nr = n;

    Value m;
    value_init(m);

    evalue EV;
    value_init(EV.d);
    value_init(EV.x.n);
    value_set_si(EV.x.n, 1);

    for (n = 0; n < nr; ++n) {
	value_assign(m, f->p[n][nc-1]);
	if (value_one_p(m) || value_mone_p(m))
	    continue;

	int j = normal_mod(f->p[n], nc-1, &m);
	if (j == nc-1) {
	    free_evalue_refs(factor);
	    value_init(factor->d);
	    evalue_set_si(factor, 0, 1);
	    break;
	}
	vec_ZZ row;
	values2zz(f->p[n], row, nc-1);
	ZZ g;
	value2zz(m, g);
	if (j < (nc-1)-1 && row[j] > g/2) {
	    for (int k = j; k < (nc-1); ++k)
		if (row[k] != 0)
		    row[k] = g - row[k];
	}

	value_init(EP.d);
	value_set_si(EP.d, 0);
	EP.x.p = new_enode(relation, 2, 0);
	value_clear(EP.x.p->arr[1].d);
	EP.x.p->arr[1] = *factor;
	evalue *ev = &EP.x.p->arr[0];
	value_set_si(ev->d, 0);
	ev->x.p = new_enode(fractional, 3, -1);
	evalue_set_si(&ev->x.p->arr[1], 0, 1);
	evalue_set_si(&ev->x.p->arr[2], 1, 1);
	evalue *E = multi_monom(row);
	value_assign(EV.d, m);
	emul(&EV, E);
	value_clear(ev->x.p->arr[0].d);
	ev->x.p->arr[0] = *E;
	delete E;
	*factor = EP;
    }

    value_clear(m);
    free_evalue_refs(&EV); 
}
#else
/*
 * 
 */
static void mask(Matrix *f, evalue *factor)
{
    int nr = f->NbRows, nc = f->NbColumns;
    int n;
    bool found = false;
    for (n = 0; n < nr && value_notzero_p(f->p[n][nc-1]); ++n)
	if (value_notone_p(f->p[n][nc-1]) &&
	    value_notmone_p(f->p[n][nc-1]))
		found = true;
    if (!found)
	return;

    Value tmp;
    value_init(tmp);
    nr = n;
    unsigned np = nc - 2;
    Vector *lcm = Vector_Alloc(np);
    Vector *val = Vector_Alloc(nc);
    Vector_Set(val->p, 0, nc);
    value_set_si(val->p[np], 1);
    Vector_Set(lcm->p, 1, np);
    for (n = 0; n < nr; ++n) {
	if (value_one_p(f->p[n][nc-1]) ||
	    value_mone_p(f->p[n][nc-1]))
	    continue;
	for (int j = 0; j < np; ++j)
	    if (value_notzero_p(f->p[n][j])) {
		Gcd(f->p[n][j], f->p[n][nc-1], &tmp);
		value_division(tmp, f->p[n][nc-1], tmp);
		value_lcm(tmp, lcm->p[j], &lcm->p[j]);
	    }
    }
    evalue EP;
    value_init(EP.d);
    mask_r(f, nr, lcm, 0, val, &EP);
    value_clear(tmp);
    Vector_Free(val);
    Vector_Free(lcm);
    emul(&EP,factor); 
    free_evalue_refs(&EP);
}
#endif

/* This structure encodes the power of the term in a rational generating function.
 * 
 * Either E == NULL or constant = 0
 * If E != NULL, then the power is 	    E
 * If E == NULL, then the power is 	    coeff * param[pos] + constant
 */
struct term_info {
    evalue	   *E;
    ZZ		    constant;
    ZZ		    coeff;
    int		    pos;
};

/* Returns the power of (t+1) in the term of a rational generating function,
 * i.e., the scalar product of the actual lattice point and lambda.
 * The lattice point is the unique lattice point in the fundamental parallelepiped
 * of the unimodual cone i shifted to the parametric vertex V.
 *
 * PD is the parameter domain, which, if != NULL, may be used to simply the
 * resulting expression.
 *
 * The result is returned in term.
 */
void lattice_point(
    Param_Vertices* V, Polyhedron *i, vec_ZZ& lambda, term_info* term,
    Polyhedron *PD)
{
    unsigned nparam = V->Vertex->NbColumns - 2;
    unsigned dim = i->Dimension;
    mat_ZZ vertex;
    vertex.SetDims(V->Vertex->NbRows, nparam+1);
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

	term->E = lattice_point(i, lambda, mv, lcm, PD);
	term->constant = 0;

	Matrix_Free(mv);
	value_clear(lcm);
	value_clear(tmp);
	return;
    }
    for (int i = 0; i < V->Vertex->NbRows; ++i) {
	assert(value_one_p(V->Vertex->p[i][nparam+1]));  // for now
	values2zz(V->Vertex->p[i], vertex[i], nparam+1);
    }

    vec_ZZ num;
    num = lambda * vertex;

    int p = -1;
    int nn = 0;
    for (int j = 0; j < nparam; ++j)
	if (num[j] != 0) {
	    ++nn;
	    p = j;
	}
    if (nn >= 2) {
	term->E = multi_monom(num);
	term->constant = 0;
    } else {
	term->E = NULL;
	term->constant = num[nparam];
	term->pos = p;
	if (p != -1)
	    term->coeff = num[p];
    }

    value_clear(lcm);
    value_clear(tmp);
}


struct counter : public np_base {
    vec_ZZ lambda;
    mat_ZZ rays;
    vec_ZZ vertex;
    vec_ZZ den;
    ZZ sign;
    ZZ num;
    int j;
    mpq_t count;

    counter(unsigned dim) : np_base(dim) {
	rays.SetDims(dim, dim);
	den.SetLength(dim);
	mpq_init(count);
    }

    void start(Polyhedron *P, unsigned MaxRays);

    ~counter() {
	mpq_clear(count);
    }

    virtual void handle_polar(Polyhedron *C, Value *vertex, QQ c);
};

struct OrthogonalException {} Orthogonal;

void counter::handle_polar(Polyhedron *C, Value *V, QQ c)
{
    int r = 0;
    add_rays(rays, C, &r);
    for (int k = 0; k < dim; ++k) {
	if (lambda * rays[k] == 0)
	    throw Orthogonal;
    }

    assert(c.d == 1);
    assert(c.n == 1 || c.n == -1);
    sign = c.n;

    lattice_point(V, C, vertex);
    num = vertex * lambda;
    den = rays * lambda;
    normalize(sign, num, den);

    dpoly d(dim, num);
    dpoly n(dim, den[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den[k], 1);
	n *= fact;
    }
    d.div(n, count, sign);
}

void counter::start(Polyhedron *P, unsigned MaxRays)
{
    for (;;) {
	try {
	    randomvector(P, lambda, dim);
	    np_base::start(P, MaxRays);
	    break;
	} catch (OrthogonalException &e) {
	    mpq_set_si(count, 0, 0);
	}
    }
}

struct bfe_term : public bfc_term_base {
    vector<evalue *> factors;

    bfe_term(int len) : bfc_term_base(len) {
    }

    ~bfe_term() {
	for (int i = 0; i < factors.size(); ++i) {
	    if (!factors[i])
		continue;
	    free_evalue_refs(factors[i]);
	    delete factors[i];
	}
    }
};

static void print_int_vector(int *v, int len, char *name)
{
    cerr << name << endl;
    for (int j = 0; j < len; ++j) {
	cerr << v[j] << " ";
    }
    cerr << endl;
}

static void print_bfc_terms(mat_ZZ& factors, bfc_vec& v)
{
    cerr << endl;
    cerr << "factors" << endl;
    cerr << factors << endl;
    for (int i = 0; i < v.size(); ++i) {
	cerr << "term: " << i << endl;
	print_int_vector(v[i]->powers, factors.NumRows(), "powers");
	cerr << "terms" << endl;
	cerr << v[i]->terms << endl;
	bfc_term* bfct = static_cast<bfc_term *>(v[i]);
	cerr << bfct->c << endl;
    }
}

static void print_bfe_terms(mat_ZZ& factors, bfc_vec& v)
{
    cerr << endl;
    cerr << "factors" << endl;
    cerr << factors << endl;
    for (int i = 0; i < v.size(); ++i) {
	cerr << "term: " << i << endl;
	print_int_vector(v[i]->powers, factors.NumRows(), "powers");
	cerr << "terms" << endl;
	cerr << v[i]->terms << endl;
	bfe_term* bfet = static_cast<bfe_term *>(v[i]);
	for (int j = 0; j < v[i]->terms.NumRows(); ++j) {
           char * test[] = {"a", "b"};
           print_evalue(stderr, bfet->factors[j], test);
           fprintf(stderr, "\n");
	}
    }
}

struct bfcounter : public bfcounter_base {
    mpq_t count;

    bfcounter(unsigned dim) : bfcounter_base(dim) {
	mpq_init(count);
	lower = 1;
    }
    ~bfcounter() {
	mpq_clear(count);
    }
    virtual void base(mat_ZZ& factors, bfc_vec& v);
};

void bfcounter::base(mat_ZZ& factors, bfc_vec& v)
{
    unsigned nf = factors.NumRows();

    for (int i = 0; i < v.size(); ++i) {
	bfc_term* bfct = static_cast<bfc_term *>(v[i]);
	int total_power = 0;
	// factor is always positive, so we always
	// change signs
	for (int k = 0; k < nf; ++k)
	    total_power += v[i]->powers[k];

	int j;
	for (j = 0; j < nf; ++j)
	    if (v[i]->powers[j] > 0)
		break;

	dpoly D(total_power, factors[j][0], 1);
	for (int k = 1; k < v[i]->powers[j]; ++k) {
	    dpoly fact(total_power, factors[j][0], 1);
	    D *= fact;
	}
	for ( ; ++j < nf; )
	    for (int k = 0; k < v[i]->powers[j]; ++k) {
		dpoly fact(total_power, factors[j][0], 1);
		D *= fact;
	    }

	for (int k = 0; k < v[i]->terms.NumRows(); ++k) {
	    dpoly n(total_power, v[i]->terms[k][0]);
	    mpq_set_si(tcount, 0, 1);
	    n.div(D, tcount, one);
	    if (total_power % 2)
		bfct->c[k].n = -bfct->c[k].n;
	    zz2value(bfct->c[k].n, tn);
	    zz2value(bfct->c[k].d, td);

	    mpz_mul(mpq_numref(tcount), mpq_numref(tcount), tn);
	    mpz_mul(mpq_denref(tcount), mpq_denref(tcount), td);
	    mpq_canonicalize(tcount);
	    mpq_add(count, count, tcount);
	}
	delete v[i];
    }
}


/* Check whether the polyhedron is unbounded and if so,
 * check whether it has any (and therefore an infinite number of)
 * integer points.
 * If one of the vertices is integer, then we are done.
 * Otherwise, transform the polyhedron such that one of the rays
 * is the first unit vector and cut it off at a height that ensures
 * that if the whole polyhedron has any points, then the remaining part
 * has integer points.  In particular we add the largest coefficient
 * of a ray to the highest vertex (rounded up).
 */
static bool Polyhedron_is_infinite(Polyhedron *P, Value* result, unsigned MaxRays)
{
    int r = 0;
    Matrix *M, *M2;
    Value c, tmp;
    Value g;
    bool first;
    Vector *v;
    Value offset, size;
    Polyhedron *R;

    if (P->NbBid == 0)
	for (; r < P->NbRays; ++r)
	    if (value_zero_p(P->Ray[r][P->Dimension+1]))
		break;
    if (P->NbBid == 0 && r == P->NbRays)
	return false;

#ifdef HAVE_LIBGLPK
    Vector *sample;

    sample = Polyhedron_Sample(P, MaxRays);
    if (!sample)
	value_set_si(*result, 0);
    else {
	value_set_si(*result, -1);
	Vector_Free(sample);
    }
    return true;
#endif

    for (int i = 0; i < P->NbRays; ++i)
	if (value_one_p(P->Ray[i][1+P->Dimension])) {
	    value_set_si(*result, -1);
	    return true;
	}

    value_init(g);
    v = Vector_Alloc(P->Dimension+1);
    Vector_Gcd(P->Ray[r]+1, P->Dimension, &g);
    Vector_AntiScale(P->Ray[r]+1, v->p, g, P->Dimension+1);
    M = unimodular_complete(v);
    value_set_si(M->p[P->Dimension][P->Dimension], 1);
    M2 = Transpose(M);
    Matrix_Free(M);
    P = Polyhedron_Preimage(P, M2, 0);
    Matrix_Free(M2);
    value_clear(g);
    Vector_Free(v);

    first = true;
    value_init(offset);
    value_init(size);
    value_init(tmp);
    value_set_si(size, 0);

    for (int i = 0; i < P->NbBid; ++i) {
	value_absolute(tmp, P->Ray[i][1]);
	if (value_gt(tmp, size))
	    value_assign(size, tmp);
    }
    for (int i = P->NbBid; i < P->NbRays; ++i) {
	if (value_zero_p(P->Ray[i][P->Dimension+1])) {
	    if (value_gt(P->Ray[i][1], size))
		value_assign(size, P->Ray[i][1]);
	    continue;
	}
	mpz_cdiv_q(tmp, P->Ray[i][1], P->Ray[i][P->Dimension+1]);
	if (first || value_gt(tmp, offset)) {
	    value_assign(offset, tmp);
	    first = false;
	}
    }
    value_addto(offset, offset, size);
    value_clear(size);
    value_clear(tmp);

    v = Vector_Alloc(P->Dimension+2);
    value_set_si(v->p[0], 1);
    value_set_si(v->p[1], -1);
    value_assign(v->p[1+P->Dimension], offset);
    R = AddConstraints(v->p, 1, P, MaxRays);
    Polyhedron_Free(P);
    P = R;

    value_clear(offset);
    Vector_Free(v);

    value_init(c);
    barvinok_count(P, &c, MaxRays);
    Polyhedron_Free(P);
    if (value_zero_p(c))
	value_set_si(*result, 0);
    else
	value_set_si(*result, -1);
    value_clear(c);

    return true;
}

typedef Polyhedron * Polyhedron_p;

static void barvinok_count_f(Polyhedron *P, Value* result, unsigned NbMaxCons);

void barvinok_count(Polyhedron *P, Value* result, unsigned NbMaxCons)
{
    unsigned dim;
    int allocated = 0;
    Polyhedron *Q;
    bool infinite = false;

    if (emptyQ2(P)) {
	value_set_si(*result, 0);
	return;
    }
    if (P->NbEq != 0) {
	do {
	    P = remove_equalities(P);
	    P = DomainConstraintSimplify(P, NbMaxCons);
	} while (!emptyQ(P) && P->NbEq != 0);
	if (emptyQ(P)) {
	    Polyhedron_Free(P);
	    value_set_si(*result, 0);
	    return;
	}
	allocated = 1;
    }
    if (Polyhedron_is_infinite(P, result, NbMaxCons)) {
	if (allocated)
	    Polyhedron_Free(P);
	return;
    }
    if (P->Dimension == 0) {
	/* Test whether the constraints are satisfied */
	POL_ENSURE_VERTICES(P);
	value_set_si(*result, !emptyQ(P));
	if (allocated)
	    Polyhedron_Free(P);
	return;
    }
    Q = Polyhedron_Factor(P, 0, NbMaxCons);
    if (Q) {
	if (allocated)
	    Polyhedron_Free(P);
	P = Q;
	allocated = 1;
    }

    barvinok_count_f(P, result, NbMaxCons);
    if (value_neg_p(*result))
	infinite = true;
    if (Q && P->next && value_notzero_p(*result)) {
	Value factor;
	value_init(factor);

	for (Q = P->next; Q; Q = Q->next) {
	    barvinok_count_f(Q, &factor, NbMaxCons);
	    if (value_neg_p(factor)) {
		infinite = true;
		continue;
	    } else if (Q->next && value_zero_p(factor)) {
		value_set_si(*result, 0);
		break;
	    }
	    value_multiply(*result, *result, factor);
	}

	value_clear(factor);
    }

    if (allocated)
	Domain_Free(P);
    if (infinite)
	value_set_si(*result, -1);
}

static void barvinok_count_f(Polyhedron *P, Value* result, unsigned NbMaxCons)
{
    if (emptyQ2(P)) {
	value_set_si(*result, 0);
	return;
    }

    if (P->Dimension == 1)
	return Line_Length(P, result);

    int c = P->NbConstraints;
    POL_ENSURE_FACETS(P);
    if (c != P->NbConstraints || P->NbEq != 0)
	return barvinok_count(P, result, NbMaxCons);

    POL_ENSURE_VERTICES(P);

    if (Polyhedron_is_infinite(P, result, NbMaxCons))
	return;

#ifdef USE_INCREMENTAL_BF
    bfcounter cnt(P->Dimension);
#elif defined USE_INCREMENTAL_DF
    icounter cnt(P->Dimension);
#else
    counter cnt(P->Dimension);
#endif
    cnt.start(P, NbMaxCons);

    assert(value_one_p(&cnt.count[0]._mp_den));
    value_assign(*result, &cnt.count[0]._mp_num);
}

static void uni_polynom(int param, Vector *c, evalue *EP)
{ 
    unsigned dim = c->Size-2;
    value_init(EP->d);
    value_set_si(EP->d,0);
    EP->x.p = new_enode(polynomial, dim+1, param+1);
    for (int j = 0; j <= dim; ++j)
	evalue_set(&EP->x.p->arr[j], c->p[j], c->p[dim+1]);
}

static void multi_polynom(Vector *c, evalue* X, evalue *EP)
{
    unsigned dim = c->Size-2;
    evalue EC;

    value_init(EC.d);
    evalue_set(&EC, c->p[dim], c->p[dim+1]);

    value_init(EP->d);
    evalue_set(EP, c->p[dim], c->p[dim+1]);
        
    for (int i = dim-1; i >= 0; --i) {
	emul(X, EP);       
	value_assign(EC.x.n, c->p[i]);
	eadd(&EC, EP);
    }
    free_evalue_refs(&EC);
}

Polyhedron *unfringe (Polyhedron *P, unsigned MaxRays)
{
    int len = P->Dimension+2;
    Polyhedron *T, *R = P;
    Value g;
    value_init(g);
    Vector *row = Vector_Alloc(len);
    value_set_si(row->p[0], 1);

    R = DomainConstraintSimplify(Polyhedron_Copy(P), MaxRays);

    Matrix *M = Matrix_Alloc(2, len-1);
    value_set_si(M->p[1][len-2], 1);
    for (int v = 0; v < P->Dimension; ++v) {
	value_set_si(M->p[0][v], 1);
	Polyhedron *I = Polyhedron_Image(R, M, 2+1);
	value_set_si(M->p[0][v], 0);
	for (int r = 0; r < I->NbConstraints; ++r) {
	    if (value_zero_p(I->Constraint[r][0]))
		continue;
	    if (value_zero_p(I->Constraint[r][1]))
		continue;
	    if (value_one_p(I->Constraint[r][1]))
		continue;
	    if (value_mone_p(I->Constraint[r][1]))
		continue;
	    value_absolute(g, I->Constraint[r][1]);
	    Vector_Set(row->p+1, 0, len-2);
	    value_division(row->p[1+v], I->Constraint[r][1], g);
	    mpz_fdiv_q(row->p[len-1], I->Constraint[r][2], g);
	    T = R;
	    R = AddConstraints(row->p, 1, R, MaxRays);
	    if (T != P)
		Polyhedron_Free(T);
	}
	Polyhedron_Free(I);
    }
    Matrix_Free(M);
    Vector_Free(row);
    value_clear(g);
    return R;
}

/* this procedure may have false negatives */
static bool Polyhedron_is_infinite_param(Polyhedron *P, unsigned nparam)
{
    int r;
    for (r = 0; r < P->NbRays; ++r) {
	if (!value_zero_p(P->Ray[r][0]) &&
		!value_zero_p(P->Ray[r][P->Dimension+1]))
	    continue;
	if (First_Non_Zero(P->Ray[r]+1+P->Dimension-nparam, nparam) == -1)
	    return true;
    }
    return false;
}

/* Check whether all rays point in the positive directions
 * for the parameters
 */
static bool Polyhedron_has_positive_rays(Polyhedron *P, unsigned nparam)
{
    int r;
    for (r = 0; r < P->NbRays; ++r)
	if (value_zero_p(P->Ray[r][P->Dimension+1])) {
	    int i;
	    for (i = P->Dimension - nparam; i < P->Dimension; ++i)
		if (value_neg_p(P->Ray[r][i+1]))
		    return false;
	}
    return true;
}

typedef evalue * evalue_p;

struct enumerator : public polar_decomposer {
    vec_ZZ lambda;
    unsigned dim, nbV;
    evalue ** vE;
    int _i;
    mat_ZZ rays;
    vec_ZZ den;
    ZZ sign;
    Polyhedron *P;
    Param_Vertices *V;
    term_info num;
    Vector *c;
    mpq_t count;

    enumerator(Polyhedron *P, unsigned dim, unsigned nbV) {
	this->P = P;
	this->dim = dim;
	this->nbV = nbV;
	randomvector(P, lambda, dim);
	rays.SetDims(dim, dim);
	den.SetLength(dim);
	c = Vector_Alloc(dim+2);

	vE = new evalue_p[nbV];
	for (int j = 0; j < nbV; ++j)
	    vE[j] = 0;

	mpq_init(count);
    }

    void decompose_at(Param_Vertices *V, int _i, unsigned MaxRays) {
	Polyhedron *C = supporting_cone_p(P, V);
	this->_i = _i;
	this->V = V;

	vE[_i] = new evalue;
	value_init(vE[_i]->d);
	evalue_set_si(vE[_i], 0, 1);

	decompose(C, MaxRays);
    }

    ~enumerator() {
	mpq_clear(count);
	Vector_Free(c);

	for (int j = 0; j < nbV; ++j)
	    if (vE[j]) {
		free_evalue_refs(vE[j]);
		delete vE[j];
	    }
	delete [] vE;
    }

    virtual void handle_polar(Polyhedron *P, int sign);
};

void enumerator::handle_polar(Polyhedron *C, int s)
{
    int r = 0;
    assert(C->NbRays-1 == dim);
    add_rays(rays, C, &r);
    for (int k = 0; k < dim; ++k) {
	if (lambda * rays[k] == 0)
	    throw Orthogonal;
    }

    sign = s;

    lattice_point(V, C, lambda, &num, 0);
    den = rays * lambda;
    normalize(sign, num.constant, den);

    dpoly n(dim, den[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den[k], 1);
	n *= fact;
    }
    if (num.E != NULL) {
	ZZ one(INIT_VAL, 1);
	dpoly_n d(dim, num.constant, one);
	d.div(n, c, sign);
	evalue EV; 
	multi_polynom(c, num.E, &EV);
	eadd(&EV , vE[_i]);
	free_evalue_refs(&EV);
	free_evalue_refs(num.E);
	delete num.E; 
    } else if (num.pos != -1) {
	dpoly_n d(dim, num.constant, num.coeff);
	d.div(n, c, sign);
	evalue EV;
	uni_polynom(num.pos, c, &EV);
	eadd(&EV , vE[_i]);
	free_evalue_refs(&EV);
    } else {
	mpq_set_si(count, 0, 1);
	dpoly d(dim, num.constant);
	d.div(n, count, sign);
	evalue EV;
	value_init(EV.d);
	evalue_set(&EV, &count[0]._mp_num, &count[0]._mp_den);
	eadd(&EV , vE[_i]);
	free_evalue_refs(&EV);
    } 
}

struct enumerator_base {
    unsigned dim;
    evalue ** vE;
    evalue ** E_vertex;
    evalue mone;
    vertex_decomposer *vpd;

    enumerator_base(unsigned dim, vertex_decomposer *vpd)
    {
	this->dim = dim;
	this->vpd = vpd;

	vE = new evalue_p[vpd->nbV];
	for (int j = 0; j < vpd->nbV; ++j)
	    vE[j] = 0;

	E_vertex = new evalue_p[dim];

	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
    }

    void decompose_at(Param_Vertices *V, int _i, unsigned MaxRays/*, Polyhedron *pVD*/) {
	//this->pVD = pVD;

	vE[_i] = new evalue;
	value_init(vE[_i]->d);
	evalue_set_si(vE[_i], 0, 1);

	vpd->decompose_at_vertex(V, _i, MaxRays);
    }

    ~enumerator_base() {
    	for (int j = 0; j < vpd->nbV; ++j)
	    if (vE[j]) {
		free_evalue_refs(vE[j]);
		delete vE[j];
	    }
	delete [] vE;

	delete [] E_vertex;

	free_evalue_refs(&mone);
    }

    evalue *E_num(int i, int d) {
	return E_vertex[i + (dim-d)];
    }
};

struct cumulator {
    evalue *factor;
    evalue *v;
    dpoly_r *r;

    cumulator(evalue *factor, evalue *v, dpoly_r *r) : 
	factor(factor), v(v), r(r) {}

    void cumulate();

    virtual void add_term(int *powers, int len, evalue *f2) = 0;
};

void cumulator::cumulate()
{
    evalue cum;  // factor * 1 * E_num[0]/1 * (E_num[0]-1)/2 *...
    evalue f;
    evalue t;	// E_num[0] - (m-1)
#ifdef USE_MODULO
    evalue *cst;
#else
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
#endif

    value_init(cum.d);
    evalue_copy(&cum, factor);
    value_init(f.d);
    value_init(f.x.n);
    value_set_si(f.d, 1);
    value_set_si(f.x.n, 1);
    value_init(t.d);
    evalue_copy(&t, v);

#ifdef USE_MODULO
    for (cst = &t; value_zero_p(cst->d); ) {
	if (cst->x.p->type == fractional)
	    cst = &cst->x.p->arr[1];
	else
	    cst = &cst->x.p->arr[0];
    }
#endif

    for (int m = 0; m < r->len; ++m) {
	if (m > 0) {
	    if (m > 1) {
		value_set_si(f.d, m);
		emul(&f, &cum);
#ifdef USE_MODULO
		value_subtract(cst->x.n, cst->x.n, cst->d);
#else
		eadd(&mone, &t);
#endif
	    }
	    emul(&t, &cum);
	}
	vector< dpoly_r_term * >& current = r->c[r->len-1-m];
	for (int j = 0; j < current.size(); ++j) {
	    if (current[j]->coeff == 0)
		continue;
	    evalue *f2 = new evalue;
	    value_init(f2->d);
	    value_init(f2->x.n);
	    zz2value(current[j]->coeff, f2->x.n);
	    zz2value(r->denom, f2->d);
	    emul(&cum, f2);

	    add_term(current[j]->powers, r->dim, f2);
	}
    }
    free_evalue_refs(&f);
    free_evalue_refs(&t);
    free_evalue_refs(&cum);
#ifndef USE_MODULO
    free_evalue_refs(&mone);
#endif
}

struct E_poly_term {
    int	    *powers;
    evalue  *E;
};

struct ie_cum : public cumulator {
    vector<E_poly_term *> terms;

    ie_cum(evalue *factor, evalue *v, dpoly_r *r) : cumulator(factor, v, r) {}

    virtual void add_term(int *powers, int len, evalue *f2);
};

void ie_cum::add_term(int *powers, int len, evalue *f2)
{
    int k;
    for (k = 0; k < terms.size(); ++k) {
	if (memcmp(terms[k]->powers, powers, len * sizeof(int)) == 0) {
	    eadd(f2, terms[k]->E);
	    free_evalue_refs(f2); 
	    delete f2;
	    break;
	}
    }
    if (k >= terms.size()) {
	E_poly_term *ET = new E_poly_term;
	ET->powers = new int[len];
	memcpy(ET->powers, powers, len * sizeof(int));
	ET->E = f2;
	terms.push_back(ET);
    }
}

struct ienumerator : public polar_decomposer, public vertex_decomposer, 
		     public enumerator_base {
    //Polyhedron *pVD;
    mat_ZZ den;
    vec_ZZ vertex;
    mpq_t tcount;

    ienumerator(Polyhedron *P, unsigned dim, unsigned nbV) :
		vertex_decomposer(P, nbV, this), enumerator_base(dim, this) {
	vertex.SetLength(dim);

	den.SetDims(dim, dim);
	mpq_init(tcount);
    }

    ~ienumerator() {
	mpq_clear(tcount);
    }

    virtual void handle_polar(Polyhedron *P, int sign);
    void reduce(evalue *factor, vec_ZZ& num, mat_ZZ& den_f);
};

void ienumerator::reduce(
	evalue *factor, vec_ZZ& num, mat_ZZ& den_f)
{
    unsigned len = den_f.NumRows();  // number of factors in den
    unsigned dim = num.length();

    if (dim == 0) {
	eadd(factor, vE[vert]);
	return;
    }

    vec_ZZ den_s;
    den_s.SetLength(len);
    mat_ZZ den_r;
    den_r.SetDims(len, dim-1);

    int r, k;

    for (r = 0; r < len; ++r) {
	den_s[r] = den_f[r][0];
	for (k = 0; k <= dim-1; ++k)
	    if (k != 0)
		den_r[r][k-(k>0)] = den_f[r][k];
    }

    ZZ num_s = num[0];
    vec_ZZ num_p;
    num_p.SetLength(dim-1);
    for (k = 0 ; k <= dim-1; ++k)
	if (k != 0)
	    num_p[k-(k>0)] = num[k];

    vec_ZZ den_p;
    den_p.SetLength(len);

    ZZ one;
    one = 1;
    normalize(one, num_s, num_p, den_s, den_p, den_r);
    if (one != 1)
	emul(&mone, factor);

    int only_param = 0;
    int no_param = 0;
    for (int k = 0; k < len; ++k) {
	if (den_p[k] == 0)
	    ++no_param;
	else if (den_s[k] == 0)
	    ++only_param;
    }
    if (no_param == 0) {
	reduce(factor, num_p, den_r);
    } else {
	int k, l;
	mat_ZZ pden;
	pden.SetDims(only_param, dim-1);

	for (k = 0, l = 0; k < len; ++k)
	    if (den_s[k] == 0)
		pden[l++] = den_r[k];

	for (k = 0; k < len; ++k)
	    if (den_p[k] == 0)
		break;

	dpoly n(no_param, num_s);
	dpoly D(no_param, den_s[k], 1);
	for ( ; ++k < len; )
	    if (den_p[k] == 0) {
		dpoly fact(no_param, den_s[k], 1);
		D *= fact;
	    }

	dpoly_r * r = 0;
	// if no_param + only_param == len then all powers
	// below will be all zero
	if (no_param + only_param == len) {
	    if (E_num(0, dim) != 0)
		r = new dpoly_r(n, len);
	    else {
		mpq_set_si(tcount, 0, 1);
		one = 1;
		n.div(D, tcount, one);

		if (value_notzero_p(mpq_numref(tcount))) {
		    evalue f;
		    value_init(f.d);
		    value_init(f.x.n);
		    value_assign(f.x.n, mpq_numref(tcount));
		    value_assign(f.d, mpq_denref(tcount));
		    emul(&f, factor);
		    reduce(factor, num_p, pden);
		    free_evalue_refs(&f);
		}
		return;
	    }
	} else {
	    for (k = 0; k < len; ++k) {
		if (den_s[k] == 0 || den_p[k] == 0)
		    continue;

		dpoly pd(no_param-1, den_s[k], 1);

		int l;
		for (l = 0; l < k; ++l)
		    if (den_r[l] == den_r[k])
			break;

		if (r == 0)
		    r = new dpoly_r(n, pd, l, len);
		else {
		    dpoly_r *nr = new dpoly_r(r, pd, l, len);
		    delete r;
		    r = nr;
		}
	    }
	}
	dpoly_r *rc = r->div(D);
	delete r;
	r = rc;
	if (E_num(0, dim) == 0) {
	    int common = pden.NumRows();
	    vector< dpoly_r_term * >& final = r->c[r->len-1];
	    int rows;
	    evalue t;
	    evalue f;
	    value_init(f.d);
	    value_init(f.x.n);
	    zz2value(r->denom, f.d);
	    for (int j = 0; j < final.size(); ++j) {
		if (final[j]->coeff == 0)
		    continue;
		rows = common;
		for (int k = 0; k < r->dim; ++k) {
		    int n = final[j]->powers[k];
		    if (n == 0)
			continue;
		    pden.SetDims(rows+n, pden.NumCols());
		    for (int l = 0; l < n; ++l)
			pden[rows+l] = den_r[k];
		    rows += n;
		}
		value_init(t.d);
		evalue_copy(&t, factor);
		zz2value(final[j]->coeff, f.x.n);
		emul(&f, &t);
		reduce(&t, num_p, pden);
		free_evalue_refs(&t);
	    }
	    free_evalue_refs(&f);
	} else {
	    ie_cum cum(factor, E_num(0, dim), r);
	    cum.cumulate();

	    int common = pden.NumRows();
	    int rows;
	    for (int j = 0; j < cum.terms.size(); ++j) {
		rows = common;
		pden.SetDims(rows, pden.NumCols());
		for (int k = 0; k < r->dim; ++k) {
		    int n = cum.terms[j]->powers[k];
		    if (n == 0)
			continue;
		    pden.SetDims(rows+n, pden.NumCols());
		    for (int l = 0; l < n; ++l)
			pden[rows+l] = den_r[k];
		    rows += n;
		}
		reduce(cum.terms[j]->E, num_p, pden);
		free_evalue_refs(cum.terms[j]->E); 
		delete cum.terms[j]->E;
		delete [] cum.terms[j]->powers;
		delete cum.terms[j];
	    }
	}
	delete r;
    }
}

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 : 
	  p->type == flooring ? 1 : 0;
}

static int edegree(evalue *e)
{
    int d = 0;
    enode *p;

    if (value_notzero_p(e->d))
        return 0;

    p = e->x.p;
    int i = type_offset(p);
    if (p->size-i-1 > d)
	d = p->size - i - 1;
    for (; i < p->size; i++) {
	int d2 = edegree(&p->arr[i]);
	if (d2 > d)
	    d = d2;
    }
    return d;
}

void ienumerator::handle_polar(Polyhedron *C, int s)
{
    assert(C->NbRays-1 == dim);

    lattice_point(V, C, vertex, E_vertex);

    int r;
    for (r = 0; r < dim; ++r)
	values2zz(C->Ray[r]+1, den[r], dim);

    evalue one;
    value_init(one.d);
    evalue_set_si(&one, s, 1);
    reduce(&one, vertex, den);
    free_evalue_refs(&one);

    for (int i = 0; i < dim; ++i)
	if (E_vertex[i]) {
	    free_evalue_refs(E_vertex[i]);
	    delete E_vertex[i];
	}
}

struct bfenumerator : public vertex_decomposer, public bf_base,
		      public enumerator_base {
    evalue *factor;

    bfenumerator(Polyhedron *P, unsigned dim, unsigned nbV) : 
		    vertex_decomposer(P, nbV, this),
		    bf_base(dim), enumerator_base(dim, this) {
	lower = 0;
	factor = NULL;
    }

    ~bfenumerator() {
    }

    virtual void handle_polar(Polyhedron *P, int sign);
    virtual void base(mat_ZZ& factors, bfc_vec& v);

    bfc_term_base* new_bf_term(int len) {
	bfe_term* t = new bfe_term(len);
	return t;
    }

    virtual void set_factor(bfc_term_base *t, int k, int change) {
	bfe_term* bfet = static_cast<bfe_term *>(t);
	factor = bfet->factors[k];
	assert(factor != NULL);
	bfet->factors[k] = NULL;
	if (change)
	    emul(&mone, factor);
    }

    virtual void set_factor(bfc_term_base *t, int k, mpq_t &q, int change) {
	bfe_term* bfet = static_cast<bfe_term *>(t);
	factor = bfet->factors[k];
	assert(factor != NULL);
	bfet->factors[k] = NULL;

	evalue f;
	value_init(f.d);
	value_init(f.x.n);
	if (change)
	    value_oppose(f.x.n, mpq_numref(q));
	else
	    value_assign(f.x.n, mpq_numref(q));
	value_assign(f.d, mpq_denref(q));
	emul(&f, factor);
	free_evalue_refs(&f);
    }

    virtual void set_factor(bfc_term_base *t, int k, const QQ& c, int change) {
	bfe_term* bfet = static_cast<bfe_term *>(t);

	factor = new evalue;

	evalue f;
	value_init(f.d);
	value_init(f.x.n);
	zz2value(c.n, f.x.n);
	if (change)
	    value_oppose(f.x.n, f.x.n);
	zz2value(c.d, f.d);

	value_init(factor->d);
	evalue_copy(factor, bfet->factors[k]);
	emul(&f, factor);
	free_evalue_refs(&f);
    }

    void set_factor(evalue *f, int change) {
	if (change)
	    emul(&mone, f);
	factor = f;
    }

    virtual void insert_term(bfc_term_base *t, int i) {
	bfe_term* bfet = static_cast<bfe_term *>(t);
	int len = t->terms.NumRows()-1;	// already increased by one

	bfet->factors.resize(len+1);
	for (int j = len; j > i; --j) {
	    bfet->factors[j] = bfet->factors[j-1];
	    t->terms[j] = t->terms[j-1];
	}
	bfet->factors[i] = factor;
	factor = NULL;
    }

    virtual void update_term(bfc_term_base *t, int i) {
	bfe_term* bfet = static_cast<bfe_term *>(t);

	eadd(factor, bfet->factors[i]);
	free_evalue_refs(factor);
	delete factor;
    }

    virtual bool constant_vertex(int dim) { return E_num(0, dim) == 0; }

    virtual void cum(bf_reducer *bfr, bfc_term_base *t, int k, dpoly_r *r);
};

struct bfe_cum : public cumulator {
    bfenumerator *bfe;
    bfc_term_base *told;
    int k;
    bf_reducer *bfr;

    bfe_cum(evalue *factor, evalue *v, dpoly_r *r, bf_reducer *bfr, 
	    bfc_term_base *t, int k, bfenumerator *e) :
	    cumulator(factor, v, r), told(t), k(k),
	    bfr(bfr), bfe(e) {
    }

    virtual void add_term(int *powers, int len, evalue *f2);
};

void bfe_cum::add_term(int *powers, int len, evalue *f2)
{
    bfr->update_powers(powers, len);

    bfc_term_base * t = bfe->find_bfc_term(bfr->vn, bfr->npowers, bfr->nnf);
    bfe->set_factor(f2, bfr->l_changes % 2);
    bfe->add_term(t, told->terms[k], bfr->l_extra_num);
}

void bfenumerator::cum(bf_reducer *bfr, bfc_term_base *t, int k,
		       dpoly_r *r)
{
    bfe_term* bfet = static_cast<bfe_term *>(t);
    bfe_cum cum(bfet->factors[k], E_num(0, bfr->d), r, bfr, t, k, this);
    cum.cumulate();
}

void bfenumerator::base(mat_ZZ& factors, bfc_vec& v)
{
    for (int i = 0; i < v.size(); ++i) {
	assert(v[i]->terms.NumRows() == 1);
	evalue *factor = static_cast<bfe_term *>(v[i])->factors[0];
	eadd(factor, vE[vert]);
	delete v[i];
    }
}

void bfenumerator::handle_polar(Polyhedron *C, int s)
{
    assert(C->NbRays-1 == enumerator_base::dim);

    bfe_term* t = new bfe_term(enumerator_base::dim);
    vector< bfc_term_base * > v;
    v.push_back(t);

    t->factors.resize(1);

    t->terms.SetDims(1, enumerator_base::dim);
    lattice_point(V, C, t->terms[0], E_vertex);

    // the elements of factors are always lexpositive
    mat_ZZ   factors;
    s = setup_factors(C, factors, t, s);

    t->factors[0] = new evalue;
    value_init(t->factors[0]->d);
    evalue_set_si(t->factors[0], s, 1);
    reduce(factors, v);

    for (int i = 0; i < enumerator_base::dim; ++i)
	if (E_vertex[i]) {
	    free_evalue_refs(E_vertex[i]);
	    delete E_vertex[i];
	}
}

#ifdef HAVE_CORRECT_VERTICES
static inline Param_Polyhedron *Polyhedron2Param_SD(Polyhedron **Din,
    Polyhedron *Cin,int WS,Polyhedron **CEq,Matrix **CT)
{
    if (WS & POL_NO_DUAL)
	WS = 0;
    return Polyhedron2Param_SimplifiedDomain(Din, Cin, WS, CEq, CT);
}
#else
static Param_Polyhedron *Polyhedron2Param_SD(Polyhedron **Din,
    Polyhedron *Cin,int WS,Polyhedron **CEq,Matrix **CT)
{
    static char data[] =    "    1   0   0   0   0   1 -18 "
			    "    1   0   0 -20   0  19   1 "
			    "    1   0   1  20   0 -20  16 "
			    "    1   0   0   0   0  -1  19 "
			    "    1   0  -1   0   0   0   4 "
			    "    1   4 -20   0   0  -1  23 "
			    "    1  -4  20   0   0   1 -22 "
			    "    1   0   1   0  20 -20  16 "
			    "    1   0   0   0 -20  19   1 ";
    static int checked = 0;
    if (!checked) {
	checked = 1;
	char *p = data;
	int n, v, i;
	Matrix *M = Matrix_Alloc(9, 7);
	for (i = 0; i < 9; ++i)
	    for (int j = 0; j < 7; ++j) {
		sscanf(p, "%d%n", &v, &n);
		p += n;
		value_set_si(M->p[i][j], v);
	    }
	Polyhedron *P = Constraints2Polyhedron(M, 1024);
	Matrix_Free(M);
	Polyhedron *U = Universe_Polyhedron(1);
	Param_Polyhedron *PP = Polyhedron2Param_Domain(P, U, 1024);
	Polyhedron_Free(P);
	Polyhedron_Free(U);
	Param_Vertices *V;
	for (i = 0, V = PP->V; V; ++i, V = V->next)
	    ;
	if (PP)
	    Param_Polyhedron_Free(PP);
	if (i != 10) {
	    fprintf(stderr, "WARNING: results may be incorrect\n");
	    fprintf(stderr, 
	"WARNING: use latest version of PolyLib to remove this warning\n");
	}
    }

    return Polyhedron2Param_SimplifiedDomain(Din, Cin, WS, CEq, CT);
}
#endif

static evalue* barvinok_enumerate_ev_f(Polyhedron *P, Polyhedron* C, 
				       unsigned MaxRays);

/* Destroys C */
static evalue* barvinok_enumerate_cst(Polyhedron *P, Polyhedron* C, 
				      unsigned MaxRays)
{
    evalue *eres;

    ALLOC(evalue, eres);
    value_init(eres->d);
    value_set_si(eres->d, 0);
    eres->x.p = new_enode(partition, 2, C->Dimension);
    EVALUE_SET_DOMAIN(eres->x.p->arr[0], DomainConstraintSimplify(C, MaxRays));
    value_set_si(eres->x.p->arr[1].d, 1);
    value_init(eres->x.p->arr[1].x.n);
    if (emptyQ(P))
	value_set_si(eres->x.p->arr[1].x.n, 0);
    else
	barvinok_count(P, &eres->x.p->arr[1].x.n, MaxRays);

    return eres;
}

evalue* barvinok_enumerate_ev(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    //P = unfringe(P, MaxRays);
    Polyhedron *Corig = C;
    Polyhedron *CEq = NULL, *rVD, *CA;
    int r = 0;
    unsigned nparam = C->Dimension;
    evalue *eres;

    evalue factor;
    value_init(factor.d);
    evalue_set_si(&factor, 1, 1);

    CA = align_context(C, P->Dimension, MaxRays);
    P = DomainIntersection(P, CA, MaxRays);
    Polyhedron_Free(CA);

    /* for now */
    POL_ENSURE_FACETS(P);
    POL_ENSURE_VERTICES(P);
    POL_ENSURE_FACETS(C);
    POL_ENSURE_VERTICES(C);

    if (C->Dimension == 0 || emptyQ(P)) {
constant:
	eres = barvinok_enumerate_cst(P, CEq ? CEq : Polyhedron_Copy(C), 
				      MaxRays);
out:
	emul(&factor, eres);
	reduce_evalue(eres);
	free_evalue_refs(&factor);
	Domain_Free(P);
	if (C != Corig)
	    Polyhedron_Free(C);
	   
	return eres;
    }
    if (Polyhedron_is_infinite_param(P, nparam))
	goto constant;

    if (P->NbEq != 0) {
	Matrix *f;
	P = remove_equalities_p(P, P->Dimension-nparam, &f);
	mask(f, &factor);
	Matrix_Free(f);
    }
    if (P->Dimension == nparam) {
	CEq = P;
	P = Universe_Polyhedron(0);
	goto constant;
    }

    Polyhedron *T = Polyhedron_Factor(P, nparam, MaxRays);
    if (T || (P->Dimension == nparam+1)) {
	Polyhedron *Q;
	Polyhedron *C2;
	for (Q = T ? T : P; Q; Q = Q->next) {
	    Polyhedron *next = Q->next;
	    Q->next = NULL;

	    Polyhedron *QC = Q;
	    if (Q->Dimension != C->Dimension)
		QC = Polyhedron_Project(Q, nparam);

	    C2 = C;
	    C = DomainIntersection(C, QC, MaxRays);
	    if (C2 != Corig)
		Polyhedron_Free(C2);
	    if (QC != Q)
		Polyhedron_Free(QC);

	    Q->next = next;
	}
    }
    if (T) {
	Polyhedron_Free(P);
	P = T;
	if (T->Dimension == C->Dimension) {
	    P = T->next;
	    T->next = NULL;
	    Polyhedron_Free(T);
	}
    }

    Polyhedron *next = P->next;
    P->next = NULL;
    eres = barvinok_enumerate_ev_f(P, C, MaxRays);
    P->next = next;

    if (P->next) {
	Polyhedron *Q;
	evalue *f;

	for (Q = P->next; Q; Q = Q->next) {
	    Polyhedron *next = Q->next;
	    Q->next = NULL;

	    f = barvinok_enumerate_ev_f(Q, C, MaxRays);
	    emul(f, eres);
	    free_evalue_refs(f);
	    free(f);

	    Q->next = next;
	}
    }

    goto out;
}

static evalue* barvinok_enumerate_ev_f(Polyhedron *P, Polyhedron* C, 
				       unsigned MaxRays)
{
    unsigned nparam = C->Dimension;

    if (P->Dimension - nparam == 1)
	return ParamLine_Length(P, C, MaxRays);

    Param_Polyhedron *PP = NULL;
    Polyhedron *CEq = NULL, *pVD;
    Matrix *CT = NULL;
    Param_Domain *D, *next;
    Param_Vertices *V;
    evalue *eres;
    Polyhedron *Porig = P;

    PP = Polyhedron2Param_SD(&P,C,MaxRays,&CEq,&CT);

    if (isIdentity(CT)) {
	Matrix_Free(CT);
	CT = NULL;
    } else {
	assert(CT->NbRows != CT->NbColumns);
	if (CT->NbRows == 1) {		// no more parameters
	    eres = barvinok_enumerate_cst(P, CEq, MaxRays);
out:
	    if (CT)
		Matrix_Free(CT);
	    if (PP)
		Param_Polyhedron_Free(PP);
	    if (P != Porig)
		Polyhedron_Free(P);

	    return eres;
	}
	nparam = CT->NbRows - 1;
    }

    unsigned dim = P->Dimension - nparam;

    ALLOC(evalue, eres);
    value_init(eres->d);
    value_set_si(eres->d, 0);

    int nd;
    for (nd = 0, D=PP->D; D; ++nd, D=D->next);
    struct section { Polyhedron *D; evalue E; };
    section *s = new section[nd];
    Polyhedron **fVD = new Polyhedron_p[nd];

try_again:
#ifdef USE_INCREMENTAL_BF
    bfenumerator et(P, dim, PP->nbV);
#elif defined USE_INCREMENTAL_DF
    ienumerator et(P, dim, PP->nbV);
#else
    enumerator et(P, dim, PP->nbV);
#endif

    for(nd = 0, D=PP->D; D; D=next) {
	next = D->next;

	Polyhedron *rVD = reduce_domain(D->Domain, CT, CEq,
					fVD, nd, MaxRays);
	if (!rVD)
	    continue;

	pVD = CT ? DomainImage(rVD,CT,MaxRays) : rVD;

	value_init(s[nd].E.d);
	evalue_set_si(&s[nd].E, 0, 1);
	s[nd].D = rVD;

	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    if (!et.vE[_i])
		try {
		    et.decompose_at(V, _i, MaxRays);
		} catch (OrthogonalException &e) {
		    if (rVD != pVD)
			Domain_Free(pVD);
		    for (; nd >= 0; --nd) {
			free_evalue_refs(&s[nd].E);
			Domain_Free(s[nd].D);
			Domain_Free(fVD[nd]);
		    }
		    goto try_again;
		}
	    eadd(et.vE[_i] , &s[nd].E);
	END_FORALL_PVertex_in_ParamPolyhedron;
	evalue_range_reduction_in_domain(&s[nd].E, pVD);

	if (CT)
	    addeliminatedparams_evalue(&s[nd].E, CT);
	++nd;
	if (rVD != pVD)
	    Domain_Free(pVD);
    }

    if (nd == 0)
	evalue_set_si(eres, 0, 1);
    else {
	eres->x.p = new_enode(partition, 2*nd, C->Dimension);
	for (int j = 0; j < nd; ++j) {
	    EVALUE_SET_DOMAIN(eres->x.p->arr[2*j], s[j].D);
	    value_clear(eres->x.p->arr[2*j+1].d);
	    eres->x.p->arr[2*j+1] = s[j].E;
	    Domain_Free(fVD[j]);
	}
    }
    delete [] s;
    delete [] fVD;

    if (CEq)
	Polyhedron_Free(CEq);
    goto out;
}

Enumeration* barvinok_enumerate(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    evalue *EP = barvinok_enumerate_ev(P, C, MaxRays);

    return partition2enumeration(EP);
}

static void SwapColumns(Value **V, int n, int i, int j)
{
    for (int r = 0; r < n; ++r)
	value_swap(V[r][i], V[r][j]);
}

static void SwapColumns(Polyhedron *P, int i, int j)
{
    SwapColumns(P->Constraint, P->NbConstraints, i, j);
    SwapColumns(P->Ray, P->NbRays, i, j);
}

/* Construct a constraint c from constraints l and u such that if
 * if constraint c holds then for each value of the other variables
 * there is at most one value of variable pos (position pos+1 in the constraints).
 *
 * Given a lower and an upper bound
 *	    n_l v_i + <c_l,x> + c_l >= 0
 *	   -n_u v_i + <c_u,x> + c_u >= 0
 * the constructed constraint is
 *
 *	    -(n_l<c_u,x> + n_u<c_l,x>) + (-n_l c_u - n_u c_l + n_l n_u - 1)
 *
 * which is then simplified to remove the content of the non-constant coefficients
 *
 * len is the total length of the constraints.
 * v is a temporary variable that can be used by this procedure
 */
static void negative_test_constraint(Value *l, Value *u, Value *c, int pos,
				     int len, Value *v)
{
    value_oppose(*v, u[pos+1]);
    Vector_Combine(l+1, u+1, c+1, *v, l[pos+1], len-1);
    value_multiply(*v, *v, l[pos+1]);
    value_subtract(c[len-1], c[len-1], *v);
    value_set_si(*v, -1);
    Vector_Scale(c+1, c+1, *v, len-1);
    value_decrement(c[len-1], c[len-1]);
    ConstraintSimplify(c, c, len, v);
}

static bool parallel_constraints(Value *l, Value *u, Value *c, int pos,
				 int len)
{
    bool parallel;
    Value g1;
    Value g2;
    value_init(g1);
    value_init(g2);

    Vector_Gcd(&l[1+pos], len, &g1);
    Vector_Gcd(&u[1+pos], len, &g2);
    Vector_Combine(l+1+pos, u+1+pos, c+1, g2, g1, len);
    parallel = First_Non_Zero(c+1, len) == -1;

    value_clear(g1);
    value_clear(g2);

    return parallel;
}

static void negative_test_constraint7(Value *l, Value *u, Value *c, int pos,
				      int exist, int len, Value *v)
{
    Value g;
    value_init(g);

    Vector_Gcd(&u[1+pos], exist, v);
    Vector_Gcd(&l[1+pos], exist, &g);
    Vector_Combine(l+1, u+1, c+1, *v, g, len-1);
    value_multiply(*v, *v, g);
    value_subtract(c[len-1], c[len-1], *v);
    value_set_si(*v, -1);
    Vector_Scale(c+1, c+1, *v, len-1);
    value_decrement(c[len-1], c[len-1]);
    ConstraintSimplify(c, c, len, v);

    value_clear(g);
}

/* Turns a x + b >= 0 into a x + b <= -1
 *
 * len is the total length of the constraint.
 * v is a temporary variable that can be used by this procedure
 */
static void oppose_constraint(Value *c, int len, Value *v)
{
    value_set_si(*v, -1);
    Vector_Scale(c+1, c+1, *v, len-1);
    value_decrement(c[len-1], c[len-1]);
}

/* Split polyhedron P into two polyhedra *pos and *neg, where
 * existential variable i has at most one solution for each
 * value of the other variables in *neg.
 *
 * The splitting is performed using constraints l and u.
 *
 * nvar: number of set variables
 * row: temporary vector that can be used by this procedure
 * f: temporary value that can be used by this procedure
 */
static bool SplitOnConstraint(Polyhedron *P, int i, int l, int u,
			      int nvar, int MaxRays, Vector *row, Value& f,
			      Polyhedron **pos, Polyhedron **neg)
{
    negative_test_constraint(P->Constraint[l], P->Constraint[u],
			     row->p, nvar+i, P->Dimension+2, &f);
    *neg = AddConstraints(row->p, 1, P, MaxRays);

    /* We found an independent, but useless constraint
     * Maybe we should detect this earlier and not
     * mark the variable as INDEPENDENT
     */
    if (emptyQ((*neg))) {
	Polyhedron_Free(*neg);
	return false;
    }

    oppose_constraint(row->p, P->Dimension+2, &f);
    *pos = AddConstraints(row->p, 1, P, MaxRays);

    if (emptyQ((*pos))) {
	Polyhedron_Free(*neg);
	Polyhedron_Free(*pos);
	return false;
    }

    return true;
}

/*
 * unimodularly transform P such that constraint r is transformed
 * into a constraint that involves only a single (the first)
 * existential variable
 *
 */
static Polyhedron *rotate_along(Polyhedron *P, int r, int nvar, int exist,
				unsigned MaxRays)
{
    Value g;
    value_init(g);

    Vector *row = Vector_Alloc(exist);
    Vector_Copy(P->Constraint[r]+1+nvar, row->p, exist);
    Vector_Gcd(row->p, exist, &g);
    if (value_notone_p(g))
	Vector_AntiScale(row->p, row->p, g, exist);
    value_clear(g);

    Matrix *M = unimodular_complete(row);
    Matrix *M2 = Matrix_Alloc(P->Dimension+1, P->Dimension+1);
    for (r = 0; r < nvar; ++r)
	value_set_si(M2->p[r][r], 1);
    for ( ; r < nvar+exist; ++r)
	Vector_Copy(M->p[r-nvar], M2->p[r]+nvar, exist);
    for ( ; r < P->Dimension+1; ++r)
	value_set_si(M2->p[r][r], 1);
    Polyhedron *T = Polyhedron_Image(P, M2, MaxRays);

    Matrix_Free(M2);
    Matrix_Free(M);
    Vector_Free(row);

    return T;
}

/* Split polyhedron P into two polyhedra *pos and *neg, where
 * existential variable i has at most one solution for each
 * value of the other variables in *neg.
 *
 * If independent is set, then the two constraints on which the
 * split will be performed need to be independent of the other
 * existential variables.
 *
 * Return true if an appropriate split could be performed.
 *
 * nvar: number of set variables
 * exist: number of existential variables
 * row: temporary vector that can be used by this procedure
 * f: temporary value that can be used by this procedure
 */
static bool SplitOnVar(Polyhedron *P, int i, 
			      int nvar, int exist, int MaxRays,
			      Vector *row, Value& f, bool independent,
			      Polyhedron **pos, Polyhedron **neg)
{
    int j;

    for (int l = P->NbEq; l < P->NbConstraints; ++l) {
	if (value_negz_p(P->Constraint[l][nvar+i+1]))
	    continue;

	if (independent) {
	    for (j = 0; j < exist; ++j)
		if (j != i && value_notzero_p(P->Constraint[l][nvar+j+1]))
		    break;
	    if (j < exist)
		continue;
	}

	for (int u = P->NbEq; u < P->NbConstraints; ++u) {
	    if (value_posz_p(P->Constraint[u][nvar+i+1]))
		continue;

	    if (independent) {
		for (j = 0; j < exist; ++j)
		    if (j != i && value_notzero_p(P->Constraint[u][nvar+j+1]))
			break;
		if (j < exist)
		    continue;
	    }

	    if (SplitOnConstraint(P, i, l, u, nvar, MaxRays, row, f, pos, neg)) {
		if (independent) {
		    if (i != 0)
			SwapColumns(*neg, nvar+1, nvar+1+i);
		}
		return true;
	    }
	}
    }

    return false;
}

static bool double_bound_pair(Polyhedron *P, int nvar, int exist,
			 int i, int l1, int l2,
			 Polyhedron **pos, Polyhedron **neg)
{
    Value f;
    value_init(f);
    Vector *row = Vector_Alloc(P->Dimension+2);
    value_set_si(row->p[0], 1);
    value_oppose(f, P->Constraint[l1][nvar+i+1]);
    Vector_Combine(P->Constraint[l1]+1, P->Constraint[l2]+1,
		   row->p+1,
		   P->Constraint[l2][nvar+i+1], f,
		   P->Dimension+1);
    ConstraintSimplify(row->p, row->p, P->Dimension+2, &f);
    *pos = AddConstraints(row->p, 1, P, 0);
    value_set_si(f, -1);
    Vector_Scale(row->p+1, row->p+1, f, P->Dimension+1);
    value_decrement(row->p[P->Dimension+1], row->p[P->Dimension+1]);
    *neg = AddConstraints(row->p, 1, P, 0);
    Vector_Free(row);
    value_clear(f);

    return !emptyQ((*pos)) && !emptyQ((*neg));
}

static bool double_bound(Polyhedron *P, int nvar, int exist,
			 Polyhedron **pos, Polyhedron **neg)
{
    for (int i = 0; i < exist; ++i) {
	int l1, l2;
	for (l1 = P->NbEq; l1 < P->NbConstraints; ++l1) {
	    if (value_negz_p(P->Constraint[l1][nvar+i+1]))
		continue;
	    for (l2 = l1 + 1; l2 < P->NbConstraints; ++l2) {
		if (value_negz_p(P->Constraint[l2][nvar+i+1]))
		    continue;
		if (double_bound_pair(P, nvar, exist, i, l1, l2, pos, neg))
		    return true;
	    }
	}
	for (l1 = P->NbEq; l1 < P->NbConstraints; ++l1) {
	    if (value_posz_p(P->Constraint[l1][nvar+i+1]))
		continue;
	    if (l1 < P->NbConstraints)
		for (l2 = l1 + 1; l2 < P->NbConstraints; ++l2) {
		    if (value_posz_p(P->Constraint[l2][nvar+i+1]))
			continue;
		if (double_bound_pair(P, nvar, exist, i, l1, l2, pos, neg))
		    return true;
	    }
	}
	return false;
    }
    return false;
}

enum constraint { 
ALL_POS = 1 << 0,
ONE_NEG = 1 << 1,
INDEPENDENT = 1 << 2,
ROT_NEG = 1 << 3
};

static evalue* enumerate_or(Polyhedron *D,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Or\n");
#endif /* DEBUG_ER */

    Polyhedron *N = D->next;
    D->next = 0;
    evalue *EP = 
	barvinok_enumerate_e(D, exist, nparam, MaxRays);
    Polyhedron_Free(D);

    for (D = N; D; D = N) {
	N = D->next;
	D->next = 0;

	evalue *EN = 
	    barvinok_enumerate_e(D, exist, nparam, MaxRays);

	eor(EN, EP);
	free_evalue_refs(EN); 
	free(EN);
	Polyhedron_Free(D);
    }

    reduce_evalue(EP);

    return EP;
}

static evalue* enumerate_sum(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    int nvar = P->Dimension - exist - nparam;
    int toswap = nvar < exist ? nvar : exist;
    for (int i = 0; i < toswap; ++i)
	SwapColumns(P, 1 + i, nvar+exist - i);
    nparam += nvar;

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Sum\n");
#endif /* DEBUG_ER */

    evalue *EP = barvinok_enumerate_e(P, exist, nparam, MaxRays);

    for (int i = 0; i < /* nvar */ nparam; ++i) {
	Matrix *C = Matrix_Alloc(1, 1 + nparam + 1);
	value_set_si(C->p[0][0], 1);
	evalue split;
	value_init(split.d);
	value_set_si(split.d, 0);
	split.x.p = new_enode(partition, 4, nparam);
	value_set_si(C->p[0][1+i], 1);
	Matrix *C2 = Matrix_Copy(C);
	EVALUE_SET_DOMAIN(split.x.p->arr[0],
	    Constraints2Polyhedron(C2, MaxRays));
	Matrix_Free(C2);
	evalue_set_si(&split.x.p->arr[1], 1, 1);
	value_set_si(C->p[0][1+i], -1);
	value_set_si(C->p[0][1+nparam], -1);
	EVALUE_SET_DOMAIN(split.x.p->arr[2],
	    Constraints2Polyhedron(C, MaxRays));
	evalue_set_si(&split.x.p->arr[3], 1, 1);
	emul(&split, EP);
	free_evalue_refs(&split);
	Matrix_Free(C);
    }
    reduce_evalue(EP);
    evalue_range_reduction(EP);

    evalue_frac2floor(EP);

    evalue *sum = esum(EP, nvar);

    free_evalue_refs(EP); 
    free(EP);
    EP = sum;

    evalue_range_reduction(EP);

    return EP;
}

static evalue* split_sure(Polyhedron *P, Polyhedron *S,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    int nvar = P->Dimension - exist - nparam;

    Matrix *M = Matrix_Alloc(exist, S->Dimension+2);
    for (int i = 0; i < exist; ++i)
	value_set_si(M->p[i][nvar+i+1], 1);
    Polyhedron *O = S;
    S = DomainAddRays(S, M, MaxRays);
    Polyhedron_Free(O);
    Polyhedron *F = DomainAddRays(P, M, MaxRays);
    Polyhedron *D = DomainDifference(F, S, MaxRays);
    O = D;
    D = Disjoint_Domain(D, 0, MaxRays);
    Polyhedron_Free(F);
    Domain_Free(O);
    Matrix_Free(M);

    M = Matrix_Alloc(P->Dimension+1-exist, P->Dimension+1);
    for (int j = 0; j < nvar; ++j)
	value_set_si(M->p[j][j], 1);
    for (int j = 0; j < nparam+1; ++j)
	value_set_si(M->p[nvar+j][nvar+exist+j], 1);
    Polyhedron *T = Polyhedron_Image(S, M, MaxRays);
    evalue *EP = barvinok_enumerate_e(T, 0, nparam, MaxRays);
    Polyhedron_Free(S);
    Polyhedron_Free(T);
    Matrix_Free(M);

    for (Polyhedron *Q = D; Q; Q = Q->next) {
	Polyhedron *N = Q->next;
	Q->next = 0;
	T = DomainIntersection(P, Q, MaxRays);
	evalue *E = barvinok_enumerate_e(T, exist, nparam, MaxRays);
	eadd(E, EP);
	free_evalue_refs(E); 
	free(E);
	Polyhedron_Free(T);
	Q->next = N;
    }
    Domain_Free(D);
    return EP;
}

static evalue* enumerate_sure(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    int i;
    Polyhedron *S = P;
    int nvar = P->Dimension - exist - nparam;
    Value lcm;
    Value f;
    value_init(lcm);
    value_init(f);

    for (i = 0; i < exist; ++i) {
	Matrix *M = Matrix_Alloc(S->NbConstraints, S->Dimension+2);
	int c = 0;
	value_set_si(lcm, 1);
	for (int j = 0; j < S->NbConstraints; ++j) {
	    if (value_negz_p(S->Constraint[j][1+nvar+i]))
		continue;
	    if (value_one_p(S->Constraint[j][1+nvar+i]))
		continue;
	    value_lcm(lcm, S->Constraint[j][1+nvar+i], &lcm);
	}

	for (int j = 0; j < S->NbConstraints; ++j) {
	    if (value_negz_p(S->Constraint[j][1+nvar+i]))
		continue;
	    if (value_one_p(S->Constraint[j][1+nvar+i]))
		continue;
	    value_division(f, lcm, S->Constraint[j][1+nvar+i]);
	    Vector_Scale(S->Constraint[j], M->p[c], f, S->Dimension+2);
	    value_subtract(M->p[c][S->Dimension+1], 
			    M->p[c][S->Dimension+1],
			    lcm);
	    value_increment(M->p[c][S->Dimension+1], 
			    M->p[c][S->Dimension+1]);
	    ++c;
	}
	Polyhedron *O = S;
	S = AddConstraints(M->p[0], c, S, MaxRays);
	if (O != P)
	    Polyhedron_Free(O);
	Matrix_Free(M);
	if (emptyQ(S)) {
	    Polyhedron_Free(S);
	    value_clear(lcm);
	    value_clear(f);
	    return 0;
	}
    }
    value_clear(lcm);
    value_clear(f);

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Sure\n");
#endif /* DEBUG_ER */

    return split_sure(P, S, exist, nparam, MaxRays);
}

static evalue* enumerate_sure2(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    int nvar = P->Dimension - exist - nparam;
    int r;
    for (r = 0; r < P->NbRays; ++r)
	if (value_one_p(P->Ray[r][0]) &&
		value_one_p(P->Ray[r][P->Dimension+1]))
	    break;

    if (r >= P->NbRays)
	return 0;

    Matrix *M = Matrix_Alloc(nvar + 1 + nparam, P->Dimension+2);
    for (int i = 0; i < nvar; ++i)
	value_set_si(M->p[i][1+i], 1);
    for (int i = 0; i < nparam; ++i)
	value_set_si(M->p[i+nvar][1+nvar+exist+i], 1);
    Vector_Copy(P->Ray[r]+1+nvar, M->p[nvar+nparam]+1+nvar, exist);
    value_set_si(M->p[nvar+nparam][0], 1);
    value_set_si(M->p[nvar+nparam][P->Dimension+1], 1);
    Polyhedron * F = Rays2Polyhedron(M, MaxRays);
    Matrix_Free(M);

    Polyhedron *I = DomainIntersection(F, P, MaxRays);
    Polyhedron_Free(F);

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Sure2\n");
#endif /* DEBUG_ER */

    return split_sure(P, I, exist, nparam, MaxRays);
}

static evalue* enumerate_cyclic(Polyhedron *P,
			  unsigned exist, unsigned nparam, 
			  evalue * EP, int r, int p, unsigned MaxRays)
{
    int nvar = P->Dimension - exist - nparam;

    /* If EP in its fractional maps only contains references
     * to the remainder parameter with appropriate coefficients
     * then we could in principle avoid adding existentially
     * quantified variables to the validity domains.
     * We'd have to replace the remainder by m { p/m }
     * and multiply with an appropriate factor that is one
     * only in the appropriate range.
     * This last multiplication can be avoided if EP
     * has a single validity domain with no (further)
     * constraints on the remainder parameter
     */

    Matrix *CT = Matrix_Alloc(nparam+1, nparam+3);
    Matrix *M = Matrix_Alloc(1, 1+nparam+3);
    for (int j = 0; j < nparam; ++j)
	if (j != p)
	    value_set_si(CT->p[j][j], 1);
    value_set_si(CT->p[p][nparam+1], 1);
    value_set_si(CT->p[nparam][nparam+2], 1);
    value_set_si(M->p[0][1+p], -1);
    value_absolute(M->p[0][1+nparam], P->Ray[0][1+nvar+exist+p]);
    value_set_si(M->p[0][1+nparam+1], 1);
    Polyhedron *CEq = Constraints2Polyhedron(M, 1);
    Matrix_Free(M);
    addeliminatedparams_enum(EP, CT, CEq, MaxRays, nparam);
    Polyhedron_Free(CEq);
    Matrix_Free(CT);

    return EP;
}

static void enumerate_vd_add_ray(evalue *EP, Matrix *Rays, unsigned MaxRays)
{
    if (value_notzero_p(EP->d))
	return;

    assert(EP->x.p->type == partition);
    assert(EP->x.p->pos == EVALUE_DOMAIN(EP->x.p->arr[0])->Dimension);
    for (int i = 0; i < EP->x.p->size/2; ++i) {
	Polyhedron *D = EVALUE_DOMAIN(EP->x.p->arr[2*i]);
	Polyhedron *N = DomainAddRays(D, Rays, MaxRays);
	EVALUE_SET_DOMAIN(EP->x.p->arr[2*i], N);
	Domain_Free(D);
    }
}

static evalue* enumerate_line(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    if (P->NbBid == 0)
	return 0;

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Line\n");
#endif /* DEBUG_ER */

    int nvar = P->Dimension - exist - nparam;
    int i, j;
    for (i = 0; i < nparam; ++i)
	if (value_notzero_p(P->Ray[0][1+nvar+exist+i]))
	    break;
    assert(i < nparam);
    for (j = i+1; j < nparam; ++j)
	if (value_notzero_p(P->Ray[0][1+nvar+exist+i]))
	    break;
    assert(j >= nparam); // for now

    Matrix *M = Matrix_Alloc(2, P->Dimension+2);
    value_set_si(M->p[0][0], 1);
    value_set_si(M->p[0][1+nvar+exist+i], 1);
    value_set_si(M->p[1][0], 1);
    value_set_si(M->p[1][1+nvar+exist+i], -1);
    value_absolute(M->p[1][1+P->Dimension], P->Ray[0][1+nvar+exist+i]);
    value_decrement(M->p[1][1+P->Dimension], M->p[1][1+P->Dimension]);
    Polyhedron *S = AddConstraints(M->p[0], 2, P, MaxRays);
    evalue *EP = barvinok_enumerate_e(S, exist, nparam, MaxRays);
    Polyhedron_Free(S);
    Matrix_Free(M);

    return enumerate_cyclic(P, exist, nparam, EP, 0, i, MaxRays);
}

static int single_param_pos(Polyhedron*P, unsigned exist, unsigned nparam, 
			    int r)
{
    int nvar = P->Dimension - exist - nparam;
    if (First_Non_Zero(P->Ray[r]+1, nvar) != -1)
	return -1;
    int i = First_Non_Zero(P->Ray[r]+1+nvar+exist, nparam);
    if (i == -1)
	return -1;
    if (First_Non_Zero(P->Ray[r]+1+nvar+exist+1, nparam-i-1) != -1)
	return -1;
    return i;
}

static evalue* enumerate_remove_ray(Polyhedron *P, int r,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
#ifdef DEBUG_ER
    fprintf(stderr, "\nER: RedundantRay\n");
#endif /* DEBUG_ER */

    Value one;
    value_init(one);
    value_set_si(one, 1);
    int len = P->NbRays-1;
    Matrix *M = Matrix_Alloc(2 * len, P->Dimension+2);
    Vector_Copy(P->Ray[0], M->p[0], r * (P->Dimension+2));
    Vector_Copy(P->Ray[r+1], M->p[r], (len-r) * (P->Dimension+2));
    for (int j = 0; j < P->NbRays; ++j) {
	if (j == r)
	    continue;
	Vector_Combine(P->Ray[j], P->Ray[r], M->p[len+j-(j>r)], 
		       one, P->Ray[j][P->Dimension+1], P->Dimension+2);
    }

    P = Rays2Polyhedron(M, MaxRays);
    Matrix_Free(M);
    evalue *EP = barvinok_enumerate_e(P, exist, nparam, MaxRays);
    Polyhedron_Free(P);
    value_clear(one);

    return EP;
}

static evalue* enumerate_redundant_ray(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    assert(P->NbBid == 0);
    int nvar = P->Dimension - exist - nparam;
    Value m;
    value_init(m);

    for (int r = 0; r < P->NbRays; ++r) {
	if (value_notzero_p(P->Ray[r][P->Dimension+1]))
	    continue;
	int i1 = single_param_pos(P, exist, nparam, r);
	if (i1 == -1)
	    continue;
	for (int r2 = r+1; r2 < P->NbRays; ++r2) {
	    if (value_notzero_p(P->Ray[r2][P->Dimension+1]))
		continue;
	    int i2 = single_param_pos(P, exist, nparam, r2);
	    if (i2 == -1)
		continue;
	    if (i1 != i2)
		continue;

	    value_division(m, P->Ray[r][1+nvar+exist+i1], 
			      P->Ray[r2][1+nvar+exist+i1]);
	    value_multiply(m, m, P->Ray[r2][1+nvar+exist+i1]);
	    /* r2 divides r => r redundant */
	    if (value_eq(m, P->Ray[r][1+nvar+exist+i1])) {
		value_clear(m);
		return enumerate_remove_ray(P, r, exist, nparam, MaxRays);
	    }

	    value_division(m, P->Ray[r2][1+nvar+exist+i1], 
			      P->Ray[r][1+nvar+exist+i1]);
	    value_multiply(m, m, P->Ray[r][1+nvar+exist+i1]);
	    /* r divides r2 => r2 redundant */
	    if (value_eq(m, P->Ray[r2][1+nvar+exist+i1])) {
		value_clear(m);
		return enumerate_remove_ray(P, r2, exist, nparam, MaxRays);
	    }
	}
    }
    value_clear(m);
    return 0;
}

static Polyhedron *upper_bound(Polyhedron *P, 
                               int pos, Value *max, Polyhedron **R)
{
    Value v;
    int r;
    value_init(v);

    *R = 0;
    Polyhedron *N;
    Polyhedron *B = 0;
    for (Polyhedron *Q = P; Q; Q = N) {
	N = Q->next;
	for (r = 0; r < P->NbRays; ++r) {
	    if (value_zero_p(P->Ray[r][P->Dimension+1]) &&
		    value_pos_p(P->Ray[r][1+pos]))
		break;
	}
	if (r < P->NbRays) {
	    Q->next = *R;
	    *R = Q;
	    continue;
	} else {
	    Q->next = B;
	    B = Q;
	}
	for (r = 0; r < P->NbRays; ++r) {
	    if (value_zero_p(P->Ray[r][P->Dimension+1]))
		continue;
	    mpz_fdiv_q(v, P->Ray[r][1+pos], P->Ray[r][1+P->Dimension]);
	    if ((!Q->next && r == 0) || value_gt(v, *max))
		value_assign(*max, v);
	}
    }
    value_clear(v);
    return B;
}

static evalue* enumerate_ray(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    assert(P->NbBid == 0);
    int nvar = P->Dimension - exist - nparam;

    int r;
    for (r = 0; r < P->NbRays; ++r)
	if (value_zero_p(P->Ray[r][P->Dimension+1]))
	    break;
    if (r >= P->NbRays)
	return 0;

    int r2;
    for (r2 = r+1; r2 < P->NbRays; ++r2)
	if (value_zero_p(P->Ray[r2][P->Dimension+1]))
	    break;
    if (r2 < P->NbRays) {
	if (nvar > 0)
	    return enumerate_sum(P, exist, nparam, MaxRays);
    }

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Ray\n");
#endif /* DEBUG_ER */

    Value m;
    Value one;
    value_init(m);
    value_init(one);
    value_set_si(one, 1);
    int i = single_param_pos(P, exist, nparam, r);
    assert(i != -1); // for now;

    Matrix *M = Matrix_Alloc(P->NbRays, P->Dimension+2);
    for (int j = 0; j < P->NbRays; ++j) {
	Vector_Combine(P->Ray[j], P->Ray[r], M->p[j], 
		       one, P->Ray[j][P->Dimension+1], P->Dimension+2);
    }
    Polyhedron *S = Rays2Polyhedron(M, MaxRays);
    Matrix_Free(M);
    Polyhedron *D = DomainDifference(P, S, MaxRays);
    Polyhedron_Free(S);
    // Polyhedron_Print(stderr, P_VALUE_FMT, D);
    assert(value_pos_p(P->Ray[r][1+nvar+exist+i])); // for now
    Polyhedron *R;
    D = upper_bound(D, nvar+exist+i, &m, &R);
    assert(D);
    Domain_Free(D);

    M = Matrix_Alloc(2, P->Dimension+2);
    value_set_si(M->p[0][0], 1);
    value_set_si(M->p[1][0], 1);
    value_set_si(M->p[0][1+nvar+exist+i], -1);
    value_set_si(M->p[1][1+nvar+exist+i], 1);
    value_assign(M->p[0][1+P->Dimension], m);
    value_oppose(M->p[1][1+P->Dimension], m);
    value_addto(M->p[1][1+P->Dimension], M->p[1][1+P->Dimension], 
		P->Ray[r][1+nvar+exist+i]);
    value_decrement(M->p[1][1+P->Dimension], M->p[1][1+P->Dimension]);
    // Matrix_Print(stderr, P_VALUE_FMT, M);
    D = AddConstraints(M->p[0], 2, P, MaxRays);
    // Polyhedron_Print(stderr, P_VALUE_FMT, D);
    value_subtract(M->p[0][1+P->Dimension], M->p[0][1+P->Dimension], 
		    P->Ray[r][1+nvar+exist+i]);
    // Matrix_Print(stderr, P_VALUE_FMT, M);
    S = AddConstraints(M->p[0], 1, P, MaxRays);
    // Polyhedron_Print(stderr, P_VALUE_FMT, S);
    Matrix_Free(M);

    evalue *EP = barvinok_enumerate_e(D, exist, nparam, MaxRays);
    Polyhedron_Free(D);
    value_clear(one);
    value_clear(m);

    if (value_notone_p(P->Ray[r][1+nvar+exist+i]))
	EP = enumerate_cyclic(P, exist, nparam, EP, r, i, MaxRays);
    else {
	M = Matrix_Alloc(1, nparam+2);
	value_set_si(M->p[0][0], 1);
	value_set_si(M->p[0][1+i], 1);
	enumerate_vd_add_ray(EP, M, MaxRays);
	Matrix_Free(M);
    }

    if (!emptyQ(S)) {
	evalue *E = barvinok_enumerate_e(S, exist, nparam, MaxRays);
	eadd(E, EP);
	free_evalue_refs(E);
	free(E);
    }
    Polyhedron_Free(S);

    if (R) {
	assert(nvar == 0);
	evalue *ER = enumerate_or(R, exist, nparam, MaxRays);
	eor(ER, EP);
	free_evalue_refs(ER);
	free(ER);
    }

    return EP;
}

static evalue* enumerate_vd(Polyhedron **PA,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    Polyhedron *P = *PA;
    int nvar = P->Dimension - exist - nparam;
    Param_Polyhedron *PP = NULL;
    Polyhedron *C = Universe_Polyhedron(nparam);
    Polyhedron *CEq;
    Matrix *CT;
    Polyhedron *PR = P;
    PP = Polyhedron2Param_SimplifiedDomain(&PR,C,MaxRays,&CEq,&CT);
    Polyhedron_Free(C);

    int nd;
    Param_Domain *D, *last;
    Value c;
    value_init(c);
    for (nd = 0, D=PP->D; D; D=D->next, ++nd)
	;

    Polyhedron **VD = new Polyhedron_p[nd];
    Polyhedron **fVD = new Polyhedron_p[nd];
    for(nd = 0, D=PP->D; D; D=D->next) {
	Polyhedron *rVD = reduce_domain(D->Domain, CT, CEq,
					fVD, nd, MaxRays);
	if (!rVD)
	    continue;

	VD[nd++] = rVD;
	last = D;
    }

    evalue *EP = 0;

    if (nd == 0)
	EP = evalue_zero();

    /* This doesn't seem to have any effect */
    if (nd == 1) {
	Polyhedron *CA = align_context(VD[0], P->Dimension, MaxRays);
	Polyhedron *O = P;
	P = DomainIntersection(P, CA, MaxRays);
	if (O != *PA)
	    Polyhedron_Free(O);
	Polyhedron_Free(CA);
	if (emptyQ(P))
	    EP = evalue_zero();
    }

    if (!EP && CT->NbColumns != CT->NbRows) {
	Polyhedron *CEqr = DomainImage(CEq, CT, MaxRays);
	Polyhedron *CA = align_context(CEqr, PR->Dimension, MaxRays);
	Polyhedron *I = DomainIntersection(PR, CA, MaxRays);
	Polyhedron_Free(CEqr);
	Polyhedron_Free(CA);
#ifdef DEBUG_ER
	fprintf(stderr, "\nER: Eliminate\n");
#endif /* DEBUG_ER */
	nparam -= CT->NbColumns - CT->NbRows;
	EP = barvinok_enumerate_e(I, exist, nparam, MaxRays);
	nparam += CT->NbColumns - CT->NbRows;
	addeliminatedparams_enum(EP, CT, CEq, MaxRays, nparam);
	Polyhedron_Free(I);
    }
    if (PR != *PA)
	Polyhedron_Free(PR);
    PR = 0;

    if (!EP && nd > 1) {
#ifdef DEBUG_ER
	fprintf(stderr, "\nER: VD\n");
#endif /* DEBUG_ER */
	for (int i = 0; i < nd; ++i) {
	    Polyhedron *CA = align_context(VD[i], P->Dimension, MaxRays);
	    Polyhedron *I = DomainIntersection(P, CA, MaxRays);

	    if (i == 0)
		EP = barvinok_enumerate_e(I, exist, nparam, MaxRays);
	    else {
		evalue *E = barvinok_enumerate_e(I, exist, nparam, MaxRays);
		eadd(E, EP);
		free_evalue_refs(E); 
		free(E);
	    }
	    Polyhedron_Free(I);
	    Polyhedron_Free(CA);
	}
    }

    for (int i = 0; i < nd; ++i) {
	Polyhedron_Free(VD[i]);
	Polyhedron_Free(fVD[i]);
    }
    delete [] VD;
    delete [] fVD;
    value_clear(c);

    if (!EP && nvar == 0) {
	Value f;
	value_init(f);
	Param_Vertices *V, *V2;
	Matrix* M = Matrix_Alloc(1, P->Dimension+2);

	FORALL_PVertex_in_ParamPolyhedron(V, last, PP) {
	    bool found = false;
	    FORALL_PVertex_in_ParamPolyhedron(V2, last, PP) {
		if (V == V2) {
		    found = true;
		    continue;
		}
		if (!found)
		    continue;
		for (int i = 0; i < exist; ++i) {
		    value_oppose(f, V->Vertex->p[i][nparam+1]);
		    Vector_Combine(V->Vertex->p[i],
				   V2->Vertex->p[i],
				   M->p[0] + 1 + nvar + exist,
				   V2->Vertex->p[i][nparam+1],
				   f,
				   nparam+1);
		    int j;
		    for (j = 0; j < nparam; ++j)
			if (value_notzero_p(M->p[0][1+nvar+exist+j]))
			    break;
		    if (j >= nparam)
			continue;
		    ConstraintSimplify(M->p[0], M->p[0], 
				       P->Dimension+2, &f);
		    value_set_si(M->p[0][0], 0);
		    Polyhedron *para = AddConstraints(M->p[0], 1, P,
						      MaxRays);
		    if (emptyQ(para)) {
			Polyhedron_Free(para);
			continue;
		    }
		    Polyhedron *pos, *neg;
		    value_set_si(M->p[0][0], 1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    neg = AddConstraints(M->p[0], 1, P, MaxRays);
		    value_set_si(f, -1);
		    Vector_Scale(M->p[0]+1, M->p[0]+1, f, 
				 P->Dimension+1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, MaxRays);
		    if (emptyQ(neg) && emptyQ(pos)) {
			Polyhedron_Free(para);
			Polyhedron_Free(pos);
			Polyhedron_Free(neg);
			continue;
		    }
#ifdef DEBUG_ER
		    fprintf(stderr, "\nER: Order\n");
#endif /* DEBUG_ER */
		    EP = barvinok_enumerate_e(para, exist, nparam, MaxRays);
		    evalue *E;
		    if (!emptyQ(pos)) {
			E = barvinok_enumerate_e(pos, exist, nparam, MaxRays);
			eadd(E, EP);
			free_evalue_refs(E); 
			free(E);
		    }
		    if (!emptyQ(neg)) {
			E = barvinok_enumerate_e(neg, exist, nparam, MaxRays);
			eadd(E, EP);
			free_evalue_refs(E); 
			free(E);
		    }
		    Polyhedron_Free(para);
		    Polyhedron_Free(pos);
		    Polyhedron_Free(neg);
		    break;
		}
		if (EP)
		    break;
	    } END_FORALL_PVertex_in_ParamPolyhedron;
	    if (EP)
		break;
	} END_FORALL_PVertex_in_ParamPolyhedron;

	if (!EP) {
	    /* Search for vertex coordinate to split on */
	    /* First look for one independent of the parameters */
	    FORALL_PVertex_in_ParamPolyhedron(V, last, PP) {
		for (int i = 0; i < exist; ++i) {
		    int j;
		    for (j = 0; j < nparam; ++j)
			if (value_notzero_p(V->Vertex->p[i][j]))
			    break;
		    if (j < nparam)
			continue;
		    value_set_si(M->p[0][0], 1);
		    Vector_Set(M->p[0]+1, 0, nvar+exist);
		    Vector_Copy(V->Vertex->p[i], 
				M->p[0] + 1 + nvar + exist, nparam+1);
		    value_oppose(M->p[0][1+nvar+i], 
				 V->Vertex->p[i][nparam+1]);

		    Polyhedron *pos, *neg;
		    value_set_si(M->p[0][0], 1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    neg = AddConstraints(M->p[0], 1, P, MaxRays);
		    value_set_si(f, -1);
		    Vector_Scale(M->p[0]+1, M->p[0]+1, f, 
				 P->Dimension+1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, MaxRays);
		    if (emptyQ(neg) || emptyQ(pos)) {
			Polyhedron_Free(pos);
			Polyhedron_Free(neg);
			continue;
		    }
		    Polyhedron_Free(pos);
		    value_increment(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, MaxRays);
#ifdef DEBUG_ER
		    fprintf(stderr, "\nER: Vertex\n");
#endif /* DEBUG_ER */
		    pos->next = neg;
		    EP = enumerate_or(pos, exist, nparam, MaxRays);
		    break;
		}
		if (EP)
		    break;
	    } END_FORALL_PVertex_in_ParamPolyhedron;
	}

	if (!EP) {
	    /* Search for vertex coordinate to split on */
	    /* Now look for one that depends on the parameters */
	    FORALL_PVertex_in_ParamPolyhedron(V, last, PP) {
		for (int i = 0; i < exist; ++i) {
		    value_set_si(M->p[0][0], 1);
		    Vector_Set(M->p[0]+1, 0, nvar+exist);
		    Vector_Copy(V->Vertex->p[i], 
				M->p[0] + 1 + nvar + exist, nparam+1);
		    value_oppose(M->p[0][1+nvar+i], 
				 V->Vertex->p[i][nparam+1]);

		    Polyhedron *pos, *neg;
		    value_set_si(M->p[0][0], 1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    neg = AddConstraints(M->p[0], 1, P, MaxRays);
		    value_set_si(f, -1);
		    Vector_Scale(M->p[0]+1, M->p[0]+1, f, 
				 P->Dimension+1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, MaxRays);
		    if (emptyQ(neg) || emptyQ(pos)) {
			Polyhedron_Free(pos);
			Polyhedron_Free(neg);
			continue;
		    }
		    Polyhedron_Free(pos);
		    value_increment(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, MaxRays);
#ifdef DEBUG_ER
		    fprintf(stderr, "\nER: ParamVertex\n");
#endif /* DEBUG_ER */
		    pos->next = neg;
		    EP = enumerate_or(pos, exist, nparam, MaxRays);
		    break;
		}
		if (EP)
		    break;
	    } END_FORALL_PVertex_in_ParamPolyhedron;
	}

	Matrix_Free(M);
	value_clear(f);
    }

    if (CEq)
	Polyhedron_Free(CEq);
    if (CT)
	Matrix_Free(CT);
    if (PP)
	Param_Polyhedron_Free(PP);
    *PA = P;

    return EP;
}

#ifndef HAVE_PIPLIB
evalue *barvinok_enumerate_pip(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    return 0;
}
#else
evalue *barvinok_enumerate_pip(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    int nvar = P->Dimension - exist - nparam;
    evalue *EP = evalue_zero();
    Polyhedron *Q, *N;

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: PIP\n");
#endif /* DEBUG_ER */

    Polyhedron *D = pip_projectout(P, nvar, exist, nparam);
    for (Q = D; Q; Q = N) {
	N = Q->next;
	Q->next = 0;
	evalue *E;
	exist = Q->Dimension - nvar - nparam;
	E = barvinok_enumerate_e(Q, exist, nparam, MaxRays);
	Polyhedron_Free(Q);
	eadd(E, EP);
	free_evalue_refs(E); 
	free(E);
    }

    return EP;
}
#endif


static bool is_single(Value *row, int pos, int len)
{
    return First_Non_Zero(row, pos) == -1 && 
	   First_Non_Zero(row+pos+1, len-pos-1) == -1;
}

static evalue* barvinok_enumerate_e_r(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays);

#ifdef DEBUG_ER
static int er_level = 0;

evalue* barvinok_enumerate_e(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    fprintf(stderr, "\nER: level %i\n", er_level);

    Polyhedron_PrintConstraints(stderr, P_VALUE_FMT, P);
    fprintf(stderr, "\nE %d\nP %d\n", exist, nparam);
    ++er_level;
    P = DomainConstraintSimplify(Polyhedron_Copy(P), MaxRays);
    evalue *EP = barvinok_enumerate_e_r(P, exist, nparam, MaxRays);
    Polyhedron_Free(P);
    --er_level;
    return EP;
}
#else
evalue* barvinok_enumerate_e(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    P = DomainConstraintSimplify(Polyhedron_Copy(P), MaxRays);
    evalue *EP = barvinok_enumerate_e_r(P, exist, nparam, MaxRays);
    Polyhedron_Free(P);
    return EP;
}
#endif

static evalue* barvinok_enumerate_e_r(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    if (exist == 0) {
	Polyhedron *U = Universe_Polyhedron(nparam);
	evalue *EP = barvinok_enumerate_ev(P, U, MaxRays);
	//char *param_name[] = {"P", "Q", "R", "S", "T" };
	//print_evalue(stdout, EP, param_name);
	Polyhedron_Free(U);
	return EP;
    }

    int nvar = P->Dimension - exist - nparam;
    int len = P->Dimension + 2;

    /* for now */
    POL_ENSURE_FACETS(P);
    POL_ENSURE_VERTICES(P);

    if (emptyQ(P))
	return evalue_zero();

    if (nvar == 0 && nparam == 0) {
	evalue *EP = evalue_zero();
	barvinok_count(P, &EP->x.n, MaxRays);
	if (value_pos_p(EP->x.n))
	    value_set_si(EP->x.n, 1);
	return EP;
    }

    int r;
    for (r = 0; r < P->NbRays; ++r)
	if (value_zero_p(P->Ray[r][0]) ||
		value_zero_p(P->Ray[r][P->Dimension+1])) {
	    int i;
	    for (i = 0; i < nvar; ++i)
		if (value_notzero_p(P->Ray[r][i+1]))
		    break;
	    if (i >= nvar)
		continue;
	    for (i = nvar + exist; i < nvar + exist + nparam; ++i)
		if (value_notzero_p(P->Ray[r][i+1]))
		    break;
	    if (i >= nvar + exist + nparam)
		break;
	}
    if (r <  P->NbRays) {
	evalue *EP = evalue_zero();
	value_set_si(EP->x.n, -1);
	return EP;
    }

    int first;
    for (r = 0; r < P->NbEq; ++r)
	if ((first = First_Non_Zero(P->Constraint[r]+1+nvar, exist)) != -1)
		break;
    if (r < P->NbEq) {
	if (First_Non_Zero(P->Constraint[r]+1+nvar+first+1, 
			   exist-first-1) != -1) {
	    Polyhedron *T = rotate_along(P, r, nvar, exist, MaxRays);
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Equality\n");
#endif /* DEBUG_ER */
	    evalue *EP = barvinok_enumerate_e(T, exist-1, nparam, MaxRays);
	    Polyhedron_Free(T);
	    return EP;
	} else {
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Fixed\n");
#endif /* DEBUG_ER */
	    if (first == 0)
		return barvinok_enumerate_e(P, exist-1, nparam, MaxRays);
	    else {
		Polyhedron *T = Polyhedron_Copy(P);
		SwapColumns(T, nvar+1, nvar+1+first);
		evalue *EP = barvinok_enumerate_e(T, exist-1, nparam, MaxRays);
		Polyhedron_Free(T);
		return EP;
	    }
	}
    }

    Vector *row = Vector_Alloc(len);
    value_set_si(row->p[0], 1);

    Value f;
    value_init(f);

    enum constraint* info = new constraint[exist];
    for (int i = 0; i < exist; ++i) {
	info[i] = ALL_POS;
	for (int l = P->NbEq; l < P->NbConstraints; ++l) {
	    if (value_negz_p(P->Constraint[l][nvar+i+1]))
		continue;
	    bool l_parallel = is_single(P->Constraint[l]+nvar+1, i, exist);
	    for (int u = P->NbEq; u < P->NbConstraints; ++u) {
		if (value_posz_p(P->Constraint[u][nvar+i+1]))
		    continue;
		bool lu_parallel = l_parallel ||
			    is_single(P->Constraint[u]+nvar+1, i, exist);
		value_oppose(f, P->Constraint[u][nvar+i+1]);
		Vector_Combine(P->Constraint[l]+1, P->Constraint[u]+1, row->p+1,
			       f, P->Constraint[l][nvar+i+1], len-1);
		if (!(info[i] & INDEPENDENT)) {
		    int j;
		    for (j = 0; j < exist; ++j)
			if (j != i && value_notzero_p(row->p[nvar+j+1]))
			    break;
		    if (j == exist) {
			//printf("independent: i: %d, l: %d, u: %d\n", i, l, u);
			info[i] = (constraint)(info[i] | INDEPENDENT);
		    }
		}
		if (info[i] & ALL_POS) {
		    value_addto(row->p[len-1], row->p[len-1], 
			      P->Constraint[l][nvar+i+1]);
		    value_addto(row->p[len-1], row->p[len-1], f);
		    value_multiply(f, f, P->Constraint[l][nvar+i+1]);
		    value_subtract(row->p[len-1], row->p[len-1], f);
		    value_decrement(row->p[len-1], row->p[len-1]);
		    ConstraintSimplify(row->p, row->p, len, &f);
		    value_set_si(f, -1);
		    Vector_Scale(row->p+1, row->p+1, f, len-1);
		    value_decrement(row->p[len-1], row->p[len-1]);
		    Polyhedron *T = AddConstraints(row->p, 1, P, MaxRays);
		    if (!emptyQ(T)) {
			//printf("not all_pos: i: %d, l: %d, u: %d\n", i, l, u);
			info[i] = (constraint)(info[i] ^ ALL_POS);
		    }
		    //puts("pos remainder");
		    //Polyhedron_Print(stdout, P_VALUE_FMT, T);
		    Polyhedron_Free(T);
		}
		if (!(info[i] & ONE_NEG)) {
		    if (lu_parallel) {
			negative_test_constraint(P->Constraint[l],
						 P->Constraint[u],
						 row->p, nvar+i, len, &f);
			oppose_constraint(row->p, len, &f);
			Polyhedron *T = AddConstraints(row->p, 1, P, MaxRays);
			if (emptyQ(T)) {
			    //printf("one_neg i: %d, l: %d, u: %d\n", i, l, u);
			    info[i] = (constraint)(info[i] | ONE_NEG);
			}
			//puts("neg remainder");
			//Polyhedron_Print(stdout, P_VALUE_FMT, T);
			Polyhedron_Free(T);
		    } else if (!(info[i] & ROT_NEG)) {
			if (parallel_constraints(P->Constraint[l],
						 P->Constraint[u],
						 row->p, nvar, exist)) {
			    negative_test_constraint7(P->Constraint[l],
						     P->Constraint[u],
						     row->p, nvar, exist,
						     len, &f);
			    oppose_constraint(row->p, len, &f);
			    Polyhedron *T = AddConstraints(row->p, 1, P, MaxRays);
			    if (emptyQ(T)) {
				// printf("rot_neg i: %d, l: %d, u: %d\n", i, l, u);
				info[i] = (constraint)(info[i] | ROT_NEG);
				r = l;
			    }
			    //puts("neg remainder");
			    //Polyhedron_Print(stdout, P_VALUE_FMT, T);
			    Polyhedron_Free(T);
			}
		    }
		}
		if (!(info[i] & ALL_POS) && (info[i] & (ONE_NEG | ROT_NEG)))
		    goto next;
	    }
	}
	if (info[i] & ALL_POS)
	    break;
next:
	;
    }

    /*
    for (int i = 0; i < exist; ++i)
	printf("%i: %i\n", i, info[i]);
    */
    for (int i = 0; i < exist; ++i)
	if (info[i] & ALL_POS) {
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Positive\n");
#endif /* DEBUG_ER */
	    // Eliminate
	    // Maybe we should chew off some of the fat here
	    Matrix *M = Matrix_Alloc(P->Dimension, P->Dimension+1);
	    for (int j = 0; j < P->Dimension; ++j)
		value_set_si(M->p[j][j + (j >= i+nvar)], 1);
	    Polyhedron *T = Polyhedron_Image(P, M, MaxRays);
	    Matrix_Free(M);
	    evalue *EP = barvinok_enumerate_e(T, exist-1, nparam, MaxRays);
	    Polyhedron_Free(T);
	    value_clear(f);
	    Vector_Free(row);
	    delete [] info;
	    return EP;
	}
    for (int i = 0; i < exist; ++i)
	if (info[i] & ONE_NEG) {
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Negative\n");
#endif /* DEBUG_ER */
	    Vector_Free(row);
	    value_clear(f);
	    delete [] info;
	    if (i == 0)
		return barvinok_enumerate_e(P, exist-1, nparam, MaxRays);
	    else {
		Polyhedron *T = Polyhedron_Copy(P);
		SwapColumns(T, nvar+1, nvar+1+i);
		evalue *EP = barvinok_enumerate_e(T, exist-1, nparam, MaxRays);
		Polyhedron_Free(T);
		return EP;
	    }
	}
    for (int i = 0; i < exist; ++i)
	if (info[i] & ROT_NEG) {
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Rotate\n");
#endif /* DEBUG_ER */
	    Vector_Free(row);
	    value_clear(f);
	    delete [] info;
	    Polyhedron *T = rotate_along(P, r, nvar, exist, MaxRays);
	    evalue *EP = barvinok_enumerate_e(T, exist-1, nparam, MaxRays);
	    Polyhedron_Free(T);
	    return EP;
	}
    for (int i = 0; i < exist; ++i)
	if (info[i] & INDEPENDENT) {
	    Polyhedron *pos, *neg;

	    /* Find constraint again and split off negative part */

	    if (SplitOnVar(P, i, nvar, exist, MaxRays,
			   row, f, true, &pos, &neg)) {
#ifdef DEBUG_ER
		fprintf(stderr, "\nER: Split\n");
#endif /* DEBUG_ER */

		evalue *EP = 
		    barvinok_enumerate_e(neg, exist-1, nparam, MaxRays);
		evalue *E = 
		    barvinok_enumerate_e(pos, exist, nparam, MaxRays);
		eadd(E, EP);
		free_evalue_refs(E); 
		free(E);
		Polyhedron_Free(neg);
		Polyhedron_Free(pos);
		value_clear(f);
		Vector_Free(row);
		delete [] info;
		return EP;
	    }
	}
    delete [] info;

    Polyhedron *O = P;
    Polyhedron *F;

    evalue *EP;

    EP = enumerate_line(P, exist, nparam, MaxRays);
    if (EP)
	goto out;

    EP = barvinok_enumerate_pip(P, exist, nparam, MaxRays);
    if (EP)
	goto out;

    EP = enumerate_redundant_ray(P, exist, nparam, MaxRays);
    if (EP)
	goto out;

    EP = enumerate_sure(P, exist, nparam, MaxRays);
    if (EP)
	goto out;

    EP = enumerate_ray(P, exist, nparam, MaxRays);
    if (EP)
	goto out;

    EP = enumerate_sure2(P, exist, nparam, MaxRays);
    if (EP)
	goto out;

    F = unfringe(P, MaxRays);
    if (!PolyhedronIncludes(F, P)) {
#ifdef DEBUG_ER
	fprintf(stderr, "\nER: Fringed\n");
#endif /* DEBUG_ER */
	EP = barvinok_enumerate_e(F, exist, nparam, MaxRays);
	Polyhedron_Free(F);
	goto out;
    }
    Polyhedron_Free(F);

    if (nparam)
	EP = enumerate_vd(&P, exist, nparam, MaxRays);
    if (EP)
	goto out2;

    if (nvar != 0) {
	EP = enumerate_sum(P, exist, nparam, MaxRays);
	goto out2;
    }

    assert(nvar == 0);

    int i;
    Polyhedron *pos, *neg;
    for (i = 0; i < exist; ++i)
	if (SplitOnVar(P, i, nvar, exist, MaxRays,
		       row, f, false, &pos, &neg))
	    break;

    assert (i < exist);

    pos->next = neg;
    EP = enumerate_or(pos, exist, nparam, MaxRays);

out2:
    if (O != P)
	Polyhedron_Free(P);

out:
    value_clear(f);
    Vector_Free(row);
    return EP;
}

/*
 * remove equalities that require a "compression" of the parameters
 */
#ifndef HAVE_COMPRESS_PARMS
static Polyhedron *remove_more_equalities(Polyhedron *P, unsigned nparam,
					  Matrix **CP, unsigned MaxRays)
{
    return P;
}
#else
static Polyhedron *remove_more_equalities(Polyhedron *P, unsigned nparam,
					  Matrix **CP, unsigned MaxRays)
{
    Matrix *M, *T;
    Polyhedron *Q;

    /* compress_parms doesn't like equalities that only involve parameters */
    for (int i = 0; i < P->NbEq; ++i)
	assert(First_Non_Zero(P->Constraint[i]+1, P->Dimension-nparam) != -1);

    M = Matrix_Alloc(P->NbEq, P->Dimension+2);
    Vector_Copy(P->Constraint[0], M->p[0], P->NbEq * (P->Dimension+2));
    *CP = compress_parms(M, nparam);
    T = align_matrix(*CP, P->Dimension+1);
    Q = Polyhedron_Preimage(P, T, MaxRays);
    Polyhedron_Free(P);
    P = Q;
    P = remove_equalities_p(P, P->Dimension-nparam, NULL);
    Matrix_Free(T);
    Matrix_Free(M);
    return P;
}
#endif

gen_fun * barvinok_series(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    Matrix *CP = NULL;
    Polyhedron *CA;
    unsigned nparam = C->Dimension;
    gen_fun *gf;

    CA = align_context(C, P->Dimension, MaxRays);
    P = DomainIntersection(P, CA, MaxRays);
    Polyhedron_Free(CA);

    if (emptyQ2(P)) {
	Polyhedron_Free(P);
	return new gen_fun;
    }

    assert(!Polyhedron_is_infinite_param(P, nparam));
    assert(P->NbBid == 0);
    assert(Polyhedron_has_positive_rays(P, nparam));
    if (P->NbEq != 0)
	P = remove_equalities_p(P, P->Dimension-nparam, NULL);
    if (P->NbEq != 0)
	P = remove_more_equalities(P, nparam, &CP, MaxRays);
    assert(P->NbEq == 0);

    gf_base *red;
    red = gf_base::create(Polyhedron_Project(P, nparam), P->Dimension, nparam);
    red->start_gf(P, MaxRays);
    Polyhedron_Free(P);
    if (CP) {
	red->gf->substitute(CP);
	Matrix_Free(CP);
    }
    gf = red->gf;
    delete red;
    return gf;
}

static Polyhedron *skew_into_positive_orthant(Polyhedron *D, unsigned nparam, 
					      unsigned MaxRays)
{
    Matrix *M = NULL;
    Value tmp;
    value_init(tmp);
    for (Polyhedron *P = D; P; P = P->next) {
	POL_ENSURE_VERTICES(P);
	assert(!Polyhedron_is_infinite_param(P, nparam));
	assert(P->NbBid == 0);
	assert(Polyhedron_has_positive_rays(P, nparam));

	for (int r = 0; r < P->NbRays; ++r) {
	    if (value_notzero_p(P->Ray[r][P->Dimension+1]))
		continue;
	    for (int i = 0; i < nparam; ++i) {
		int j;
		if (value_posz_p(P->Ray[r][i+1]))
		    continue;
		if (!M) {
		    M = Matrix_Alloc(D->Dimension+1, D->Dimension+1);
		    for (int i = 0; i < D->Dimension+1; ++i)
			value_set_si(M->p[i][i], 1);
		} else {
		    Inner_Product(P->Ray[r]+1, M->p[i], D->Dimension+1, &tmp);
		    if (value_posz_p(tmp))
			continue;
		}
		for (j = P->Dimension - nparam; j < P->Dimension; ++j)
		    if (value_pos_p(P->Ray[r][j+1]))
			break;
		assert(j < P->Dimension);
		value_pdivision(tmp, P->Ray[r][j+1], P->Ray[r][i+1]);
		value_subtract(M->p[i][j], M->p[i][j], tmp);
	    }
	}
    }
    value_clear(tmp);
    if (M) {
	D = DomainImage(D, M, MaxRays);
	Matrix_Free(M);
    }
    return D;
}

gen_fun* barvinok_enumerate_union_series(Polyhedron *D, Polyhedron* C, 
					 unsigned MaxRays)
{
    Polyhedron *conv, *D2;
    gen_fun *gf = NULL, *gf2;
    unsigned nparam = C->Dimension;
    ZZ one, mone;
    one = 1;
    mone = -1;
    D2 = skew_into_positive_orthant(D, nparam, MaxRays);
    for (Polyhedron *P = D2; P; P = P->next) {
	assert(P->Dimension == D2->Dimension);
	POL_ENSURE_VERTICES(P);
	/* it doesn't matter which reducer we use, since we don't actually
	 * reduce anything here
	 */
	partial_reducer red(Polyhedron_Project(P, P->Dimension), P->Dimension, 
			    P->Dimension);
	red.start(P, MaxRays);
	if (!gf)
	    gf = red.gf;
	else {
	    gf->add_union(red.gf, MaxRays);
	    delete red.gf;
	}
    }
    /* we actually only need the convex union of the parameter space
     * but the reducer classes currently expect a polyhedron in
     * the combined space
     */
    Polyhedron_Free(gf->context);
    gf->context = DomainConvex(D2, MaxRays);

    gf2 = gf->summate(D2->Dimension - nparam);

    delete gf;
    if (D != D2)
	Domain_Free(D2);
    return gf2;
}

evalue* barvinok_enumerate_union(Polyhedron *D, Polyhedron* C, unsigned MaxRays)
{
    evalue *EP;
    gen_fun *gf = barvinok_enumerate_union_series(D, C, MaxRays);
    EP = *gf;
    delete gf;
    return EP;
}
