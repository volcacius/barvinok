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
#include <barvinok/evalue.h>
extern "C" {
#include "piputil.h"
}
#include "config.h"
#include <barvinok/barvinok.h>
#include <barvinok/genfun.h>
#include <barvinok/options.h>
#include <barvinok/sample.h>
#include "conversion.h"
#include "decomposer.h"
#include "lattice_point.h"
#include "reduce_domain.h"
#include "genfun_constructor.h"
#include "remove_equalities.h"
#include "scale.h"
#include "volume.h"

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

static void mask_fractional(Matrix *f, evalue *factor)
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

/*
 * 
 */
static void mask_table(Matrix *f, evalue *factor)
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

static void mask(Matrix *f, evalue *factor, barvinok_options *options)
{
    if (options->lookup_table)
	mask_table(f, factor);
    else
	mask_fractional(f, factor);
}

struct counter : public np_base {
    vec_ZZ lambda;
    mat_ZZ vertex;
    vec_ZZ den;
    ZZ sign;
    vec_ZZ num;
    ZZ offset;
    int j;
    mpq_t count;

    counter(unsigned dim) : np_base(dim) {
	den.SetLength(dim);
	mpq_init(count);
    }

    virtual void init(Polyhedron *P) {
	randomvector(P, lambda, dim);
    }

    virtual void reset() {
	mpq_set_si(count, 0, 0);
    }

    ~counter() {
	mpq_clear(count);
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, int *closed, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};

void counter::handle(const mat_ZZ& rays, Value *V, const QQ& c, unsigned long det,
		     int *closed, barvinok_options *options)
{
    for (int k = 0; k < dim; ++k) {
	if (lambda * rays[k] == 0)
	    throw Orthogonal;
    }

    assert(c.d == 1);
    assert(c.n == 1 || c.n == -1);
    sign = c.n;

    lattice_point(V, rays, vertex, det, closed);
    num = vertex * lambda;
    den = rays * lambda;
    offset = 0;
    normalize(sign, offset, den);

    num[0] += offset;
    dpoly d(dim, num[0]);
    for (int k = 1; k < num.length(); ++k) {
	num[k] += offset;
	dpoly term(dim, num[k]);
	d += term;
    }
    dpoly n(dim, den[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den[k], 1);
	n *= fact;
    }
    d.div(n, count, sign);
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
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
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
static bool Polyhedron_is_infinite(Polyhedron *P, Value* result,
				   barvinok_options *options)
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

    if (options->count_sample_infinite) {
	Vector *sample;

	sample = Polyhedron_Sample(P, options);
	if (!sample)
	    value_set_si(*result, 0);
	else {
	    value_set_si(*result, -1);
	    Vector_Free(sample);
	}
	return true;
    }

    for (int i = 0; i < P->NbRays; ++i)
	if (value_one_p(P->Ray[i][1+P->Dimension])) {
	    value_set_si(*result, -1);
	    return true;
	}

    value_init(g);
    M = Matrix_Alloc(P->Dimension+1, P->Dimension+1);
    Vector_Gcd(P->Ray[r]+1, P->Dimension, &g);
    Vector_AntiScale(P->Ray[r]+1, M->p[0], g, P->Dimension+1);
    int ok = unimodular_complete(M, 1);
    assert(ok);
    value_set_si(M->p[P->Dimension][P->Dimension], 1);
    M2 = Transpose(M);
    Matrix_Free(M);
    P = Polyhedron_Preimage(P, M2, 0);
    Matrix_Free(M2);
    value_clear(g);

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
    R = AddConstraints(v->p, 1, P, options->MaxRays);
    Polyhedron_Free(P);
    P = R;

    value_clear(offset);
    Vector_Free(v);

    value_init(c);
    barvinok_count_with_options(P, &c, options);
    Polyhedron_Free(P);
    if (value_zero_p(c))
	value_set_si(*result, 0);
    else
	value_set_si(*result, -1);
    value_clear(c);

    return true;
}

typedef Polyhedron * Polyhedron_p;

static void barvinok_count_f(Polyhedron *P, Value* result,
			     barvinok_options *options);

void barvinok_count_with_options(Polyhedron *P, Value* result,
				 struct barvinok_options *options)
{
    unsigned dim;
    int allocated = 0;
    Polyhedron *Q;
    bool infinite = false;

    if (P->next)
	fprintf(stderr,
	    "barvinok_count: input is a union; only first polyhedron is counted\n");

    if (emptyQ2(P)) {
	value_set_si(*result, 0);
	return;
    }
    if (P->NbEq != 0) {
	Q = NULL;
	do {
	    P = remove_equalities(P, options->MaxRays);
	    P = DomainConstraintSimplify(P, options->MaxRays);
	    if (Q)
		Polyhedron_Free(Q);
	    Q = P;
	} while (!emptyQ(P) && P->NbEq != 0);
	if (emptyQ(P)) {
	    Polyhedron_Free(P);
	    value_set_si(*result, 0);
	    return;
	}
	allocated = 1;
    }
    if (Polyhedron_is_infinite(P, result, options)) {
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
    Q = Polyhedron_Factor(P, 0, NULL, options->MaxRays);
    if (Q) {
	if (allocated)
	    Polyhedron_Free(P);
	P = Q;
	allocated = 1;
    }

    barvinok_count_f(P, result, options);
    if (value_neg_p(*result))
	infinite = true;
    if (Q && P->next && value_notzero_p(*result)) {
	Value factor;
	value_init(factor);

	for (Q = P->next; Q; Q = Q->next) {
	    barvinok_count_f(Q, &factor, options);
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

void barvinok_count(Polyhedron *P, Value* result, unsigned NbMaxCons)
{
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = NbMaxCons;
    barvinok_count_with_options(P, result, options);
    barvinok_options_free(options);
}

static void barvinok_count_f(Polyhedron *P, Value* result,
			     barvinok_options *options)
{
    if (emptyQ2(P)) {
	value_set_si(*result, 0);
	return;
    }

    if (P->Dimension == 1)
	return Line_Length(P, result);

    int c = P->NbConstraints;
    POL_ENSURE_FACETS(P);
    if (c != P->NbConstraints || P->NbEq != 0) {
	Polyhedron *next = P->next;
	P->next = NULL;
	barvinok_count_with_options(P, result, options);
	P->next = next;
	return;
    }

    POL_ENSURE_VERTICES(P);

    if (Polyhedron_is_infinite(P, result, options))
	return;

    np_base *cnt;
    if (options->incremental_specialization == 2)
	cnt = new bfcounter(P->Dimension);
    else if (options->incremental_specialization == 1)
	cnt = new icounter(P->Dimension);
    else
	cnt = new counter(P->Dimension);
    cnt->start(P, options);

    cnt->get_count(result);
    delete cnt;
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

struct enumerator_base {
    unsigned dim;
    evalue ** vE;
    evalue mone;
    vertex_decomposer *vpd;

    enumerator_base(unsigned dim, vertex_decomposer *vpd)
    {
	this->dim = dim;
	this->vpd = vpd;

	vE = new evalue_p[vpd->nbV];
	for (int j = 0; j < vpd->nbV; ++j)
	    vE[j] = 0;

	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
    }

    void decompose_at(Param_Vertices *V, int _i, barvinok_options *options) {
	//this->pVD = pVD;

	vE[_i] = new evalue;
	value_init(vE[_i]->d);
	evalue_set_si(vE[_i], 0, 1);

	vpd->decompose_at_vertex(V, _i, options);
    }

    virtual ~enumerator_base() {
    	for (int j = 0; j < vpd->nbV; ++j)
	    if (vE[j]) {
		free_evalue_refs(vE[j]);
		delete vE[j];
	    }
	delete [] vE;

	free_evalue_refs(&mone);
    }

    static enumerator_base *create(Polyhedron *P, unsigned dim, unsigned nbV,
				     barvinok_options *options);
};

struct enumerator : public signed_cone_consumer, public vertex_decomposer,
		    public enumerator_base {
    vec_ZZ lambda;
    vec_ZZ den;
    ZZ sign;
    term_info num;
    Vector *c;
    mpq_t count;

    enumerator(Polyhedron *P, unsigned dim, unsigned nbV) :
		vertex_decomposer(P, nbV, *this), enumerator_base(dim, this) {
	this->P = P;
	this->nbV = nbV;
	randomvector(P, lambda, dim);
	den.SetLength(dim);
	c = Vector_Alloc(dim+2);

	mpq_init(count);
    }

    ~enumerator() {
	mpq_clear(count);
	Vector_Free(c);
    }

    virtual void handle(const signed_cone& sc, barvinok_options *options);
};

void enumerator::handle(const signed_cone& sc, barvinok_options *options)
{
    int r = 0;
    assert(sc.rays.NumRows() == dim);
    for (int k = 0; k < dim; ++k) {
	if (lambda * sc.rays[k] == 0)
	    throw Orthogonal;
    }

    sign = sc.sign;

    lattice_point(V, sc.rays, lambda, &num, sc.det, sc.closed, options);
    den = sc.rays * lambda;
    ZZ offset;
    normalize(sign, offset, den);

    dpoly n(dim, den[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den[k], 1);
	n *= fact;
    }
    if (num.E != NULL) {
	ZZ one(INIT_VAL, 1);
	dpoly_n d(dim, offset, one);
	d.div(n, c, sign);
	for (unsigned long i = 0; i < sc.det; ++i) {
	    evalue EV;
	    multi_polynom(c, num.E[i], &EV);
	    eadd(&EV , vE[vert]);
	    free_evalue_refs(&EV);
	    free_evalue_refs(num.E[i]);
	    delete num.E[i];
	}
	delete [] num.E; 
    } else {
	mpq_set_si(count, 0, 1);
	if (num.constant.length() == 1) {
	    num.constant[0] += offset;
	    dpoly d(dim, num.constant[0]);
	    d.div(n, count, sign);
	} else {
	    ZZ one(INIT_VAL, 1);
	    dpoly_n d(dim, offset, one);
	    d.div(n, c, sign);
	    Value x, sum, acc;
	    value_init(x);
	    value_init(acc);
	    for (unsigned long i = 0; i < sc.det; ++i) {
		value_assign(acc, c->p[dim]);
		zz2value(num.constant[i], x);
		for (int j = dim-1; j >= 0; --j) {
		    value_multiply(acc, acc, x);
		    value_addto(acc, acc, c->p[j]);
		}
		value_addto(mpq_numref(count), mpq_numref(count), acc);
	    }
	    mpz_set(mpq_denref(count), c->p[dim+1]);
	    value_clear(acc);
	    value_clear(x);
	}
	evalue EV;
	value_init(EV.d);
	evalue_set(&EV, &count[0]._mp_num, &count[0]._mp_den);
	eadd(&EV, vE[vert]);
	free_evalue_refs(&EV);
    } 
}

struct ienumerator_base : enumerator_base {
    evalue ** E_vertex;

    ienumerator_base(unsigned dim, vertex_decomposer *vpd) :
			enumerator_base(dim,vpd) {
	E_vertex = new evalue_p[dim];
    }

    virtual ~ienumerator_base() {
	delete [] E_vertex;
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

    void cumulate(barvinok_options *options);

    virtual void add_term(const vector<int>& powers, evalue *f2) = 0;
    virtual ~cumulator() {}
};

void cumulator::cumulate(barvinok_options *options)
{
    evalue cum;  // factor * 1 * E_num[0]/1 * (E_num[0]-1)/2 *...
    evalue f;
    evalue t;	// E_num[0] - (m-1)
    evalue *cst;
    evalue mone;

    if (options->lookup_table) {
	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
    }

    value_init(cum.d);
    evalue_copy(&cum, factor);
    value_init(f.d);
    value_init(f.x.n);
    value_set_si(f.d, 1);
    value_set_si(f.x.n, 1);
    value_init(t.d);
    evalue_copy(&t, v);

    if (!options->lookup_table) {
	for (cst = &t; value_zero_p(cst->d); ) {
	    if (cst->x.p->type == fractional)
		cst = &cst->x.p->arr[1];
	    else
		cst = &cst->x.p->arr[0];
	}
    }

    for (int m = 0; m < r->len; ++m) {
	if (m > 0) {
	    if (m > 1) {
		value_set_si(f.d, m);
		emul(&f, &cum);
		if (!options->lookup_table)
		    value_subtract(cst->x.n, cst->x.n, cst->d);
		else
		    eadd(&mone, &t);
	    }
	    emul(&t, &cum);
	}
	dpoly_r_term_list& current = r->c[r->len-1-m];
	dpoly_r_term_list::iterator j;
	for (j = current.begin(); j != current.end(); ++j) {
	    if ((*j)->coeff == 0)
		continue;
	    evalue *f2 = new evalue;
	    value_init(f2->d);
	    value_init(f2->x.n);
	    zz2value((*j)->coeff, f2->x.n);
	    zz2value(r->denom, f2->d);
	    emul(&cum, f2);

	    add_term((*j)->powers, f2);
	}
    }
    free_evalue_refs(&f);
    free_evalue_refs(&t);
    free_evalue_refs(&cum);
    if (options->lookup_table)
	free_evalue_refs(&mone);
}

struct E_poly_term {
    vector<int>	powers;
    evalue  *E;
};

struct ie_cum : public cumulator {
    vector<E_poly_term *> terms;

    ie_cum(evalue *factor, evalue *v, dpoly_r *r) : cumulator(factor, v, r) {}

    virtual void add_term(const vector<int>& powers, evalue *f2);
};

void ie_cum::add_term(const vector<int>& powers, evalue *f2)
{
    int k;
    for (k = 0; k < terms.size(); ++k) {
	if (terms[k]->powers == powers) {
	    eadd(f2, terms[k]->E);
	    free_evalue_refs(f2); 
	    delete f2;
	    break;
	}
    }
    if (k >= terms.size()) {
	E_poly_term *ET = new E_poly_term;
	ET->powers = powers;
	ET->E = f2;
	terms.push_back(ET);
    }
}

struct ienumerator : public signed_cone_consumer, public vertex_decomposer,
		     public ienumerator_base {
    //Polyhedron *pVD;
    mat_ZZ den;
    mat_ZZ vertex;
    mpq_t tcount;

    ienumerator(Polyhedron *P, unsigned dim, unsigned nbV) :
		vertex_decomposer(P, nbV, *this), ienumerator_base(dim, this) {
	vertex.SetDims(1, dim);

	den.SetDims(dim, dim);
	mpq_init(tcount);
    }

    ~ienumerator() {
	mpq_clear(tcount);
    }

    virtual void handle(const signed_cone& sc, barvinok_options *options);
    void reduce(evalue *factor, const mat_ZZ& num, const mat_ZZ& den_f,
		barvinok_options *options);
};

void ienumerator::reduce(evalue *factor, const mat_ZZ& num, const mat_ZZ& den_f,
			 barvinok_options *options)
{
    unsigned len = den_f.NumRows();  // number of factors in den
    unsigned dim = num.NumCols();
    assert(num.NumRows() == 1);

    if (dim == 0) {
	eadd(factor, vE[vert]);
	return;
    }

    vec_ZZ den_s;
    mat_ZZ den_r;
    vec_ZZ num_s;
    mat_ZZ num_p;

    split_one(num, num_s, num_p, den_f, den_s, den_r);

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
	reduce(factor, num_p, den_r, options);
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

	dpoly n(no_param, num_s[0]);
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
		    reduce(factor, num_p, pden, options);
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
	    dpoly_r_term_list& final = r->c[r->len-1];
	    int rows;
	    evalue t;
	    evalue f;
	    value_init(f.d);
	    value_init(f.x.n);
	    zz2value(r->denom, f.d);
	    dpoly_r_term_list::iterator j;
	    for (j = final.begin(); j != final.end(); ++j) {
		if ((*j)->coeff == 0)
		    continue;
		rows = common;
		for (int k = 0; k < r->dim; ++k) {
		    int n = (*j)->powers[k];
		    if (n == 0)
			continue;
		    pden.SetDims(rows+n, pden.NumCols());
		    for (int l = 0; l < n; ++l)
			pden[rows+l] = den_r[k];
		    rows += n;
		}
		value_init(t.d);
		evalue_copy(&t, factor);
		zz2value((*j)->coeff, f.x.n);
		emul(&f, &t);
		reduce(&t, num_p, pden, options);
		free_evalue_refs(&t);
	    }
	    free_evalue_refs(&f);
	} else {
	    ie_cum cum(factor, E_num(0, dim), r);
	    cum.cumulate(options);

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
		reduce(cum.terms[j]->E, num_p, pden, options);
		free_evalue_refs(cum.terms[j]->E); 
		delete cum.terms[j]->E;
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

void ienumerator::handle(const signed_cone& sc, barvinok_options *options)
{
    assert(sc.det == 1);
    assert(!sc.closed);
    assert(sc.rays.NumRows() == dim);

    lattice_point(V, sc.rays, vertex[0], E_vertex, options);

    den = sc.rays;

    evalue one;
    value_init(one.d);
    evalue_set_si(&one, sc.sign, 1);
    reduce(&one, vertex, den, options);
    free_evalue_refs(&one);

    for (int i = 0; i < dim; ++i)
	if (E_vertex[i]) {
	    free_evalue_refs(E_vertex[i]);
	    delete E_vertex[i];
	}
}

struct bfenumerator : public vertex_decomposer, public bf_base,
		      public ienumerator_base {
    evalue *factor;

    bfenumerator(Polyhedron *P, unsigned dim, unsigned nbV) : 
		    vertex_decomposer(P, nbV, *this),
		    bf_base(dim), ienumerator_base(dim, this) {
	lower = 0;
	factor = NULL;
    }

    ~bfenumerator() {
    }

    virtual void handle(const signed_cone& sc, barvinok_options *options);
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

    virtual void cum(bf_reducer *bfr, bfc_term_base *t, int k, dpoly_r *r,
		     barvinok_options *options);
};

enumerator_base *enumerator_base::create(Polyhedron *P, unsigned dim, unsigned nbV,
					 barvinok_options *options)
{
    enumerator_base *eb;

    if (options->incremental_specialization == BV_SPECIALIZATION_BF)
	eb = new bfenumerator(P, dim, nbV);
    else if (options->incremental_specialization == BV_SPECIALIZATION_DF)
	eb = new ienumerator(P, dim, nbV);
    else
	eb = new enumerator(P, dim, nbV);

    return eb;
}

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

    virtual void add_term(const vector<int>& powers, evalue *f2);
};

void bfe_cum::add_term(const vector<int>& powers, evalue *f2)
{
    bfr->update_powers(powers);

    bfc_term_base * t = bfe->find_bfc_term(bfr->vn, bfr->npowers, bfr->nnf);
    bfe->set_factor(f2, bfr->l_changes % 2);
    bfe->add_term(t, told->terms[k], bfr->l_extra_num);
}

void bfenumerator::cum(bf_reducer *bfr, bfc_term_base *t, int k,
		       dpoly_r *r, barvinok_options *options)
{
    bfe_term* bfet = static_cast<bfe_term *>(t);
    bfe_cum cum(bfet->factors[k], E_num(0, bfr->d), r, bfr, t, k, this);
    cum.cumulate(options);
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

void bfenumerator::handle(const signed_cone& sc, barvinok_options *options)
{
    assert(sc.det == 1);
    assert(!sc.closed);
    assert(sc.rays.NumRows() == enumerator_base::dim);

    bfe_term* t = new bfe_term(enumerator_base::dim);
    vector< bfc_term_base * > v;
    v.push_back(t);

    t->factors.resize(1);

    t->terms.SetDims(1, enumerator_base::dim);
    lattice_point(V, sc.rays, t->terms[0], E_vertex, options);

    // the elements of factors are always lexpositive
    mat_ZZ   factors;
    int s = setup_factors(sc.rays, factors, t, sc.sign);

    t->factors[0] = new evalue;
    value_init(t->factors[0]->d);
    evalue_set_si(t->factors[0], s, 1);
    reduce(factors, v, options);

    for (int i = 0; i < enumerator_base::dim; ++i)
	if (E_vertex[i]) {
	    free_evalue_refs(E_vertex[i]);
	    delete E_vertex[i];
	}
}

static inline Param_Polyhedron *Polyhedron2Param_MR(Polyhedron *Din,
    Polyhedron *Cin, int WS)
{
    if (WS & POL_NO_DUAL)
	WS = 0;
    return Polyhedron2Param_Domain(Din, Cin, WS);
}

static evalue* barvinok_enumerate_ev_f(Polyhedron *P, Polyhedron* C, 
				       barvinok_options *options);

/* Destroys C */
static evalue* barvinok_enumerate_cst(Polyhedron *P, Polyhedron* C, 
				      struct barvinok_options *options)
{
    evalue *eres;

    ALLOC(evalue, eres);
    value_init(eres->d);
    value_set_si(eres->d, 0);
    eres->x.p = new_enode(partition, 2, C->Dimension);
    EVALUE_SET_DOMAIN(eres->x.p->arr[0],
		      DomainConstraintSimplify(C, options->MaxRays));
    value_set_si(eres->x.p->arr[1].d, 1);
    value_init(eres->x.p->arr[1].x.n);
    if (emptyQ2(P))
	value_set_si(eres->x.p->arr[1].x.n, 0);
    else
	barvinok_count_with_options(P, &eres->x.p->arr[1].x.n, options);

    return eres;
}

/* frees P */
static evalue* enumerate(Polyhedron *P, Polyhedron* C,
					struct barvinok_options *options)
{
    //P = unfringe(P, MaxRays);
    Polyhedron *next;
    Polyhedron *Corig = C;
    Polyhedron *CEq = NULL, *rVD;
    int r = 0;
    unsigned nparam = C->Dimension;
    evalue *eres;
    Matrix *CP = NULL;

    evalue factor;
    value_init(factor.d);
    evalue_set_si(&factor, 1, 1);

    /* for now */
    POL_ENSURE_FACETS(P);
    POL_ENSURE_VERTICES(P);
    POL_ENSURE_FACETS(C);
    POL_ENSURE_VERTICES(C);

    if (C->Dimension == 0 || emptyQ(P)) {
constant:
	eres = barvinok_enumerate_cst(P, CEq ? CEq : Polyhedron_Copy(C), options);
out:
	if (CP) {
	    evalue_backsubstitute(eres, CP, options->MaxRays);
	    Matrix_Free(CP);
	}

	emul(&factor, eres);
	if (options->approximation_method == BV_APPROX_DROP) {
	    if (options->polynomial_approximation == BV_APPROX_SIGN_UPPER)
		evalue_frac2polynomial(eres, 1, options->MaxRays);
	    if (options->polynomial_approximation == BV_APPROX_SIGN_LOWER)
		evalue_frac2polynomial(eres, -1, options->MaxRays);
	    if (options->polynomial_approximation == BV_APPROX_SIGN_APPROX)
		evalue_frac2polynomial(eres, 0, options->MaxRays);
	}
	reduce_evalue(eres);
	free_evalue_refs(&factor);
	Domain_Free(P);
	if (C != Corig)
	    Polyhedron_Free(C);
	   
	return eres;
    }
    if (Polyhedron_is_unbounded(P, nparam, options->MaxRays))
	goto constant;

    if (P->NbEq != 0) {
	Matrix *f;
	P = remove_equalities_p(P, P->Dimension-nparam, &f, options->MaxRays);
	mask(f, &factor, options);
	Matrix_Free(f);
    }
    if (P->Dimension == nparam) {
	CEq = P;
	P = Universe_Polyhedron(0);
	goto constant;
    }
    if (P->NbEq != 0) {
	Polyhedron *Q = P;
	Polyhedron *D = C;
	remove_all_equalities(&Q, &C, &CP, NULL, nparam, options->MaxRays);
	if (C != D && D != Corig)
	    Polyhedron_Free(D);
	eres = enumerate(Q, C, options);
	goto out;
    }

    Polyhedron *T = Polyhedron_Factor(P, nparam, NULL, options->MaxRays);
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
	    C = DomainIntersection(C, QC, options->MaxRays);
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

    next = P->next;
    P->next = NULL;
    eres = barvinok_enumerate_ev_f(P, C, options);
    P->next = next;

    if (P->next) {
	Polyhedron *Q;
	evalue *f;

	for (Q = P->next; Q; Q = Q->next) {
	    Polyhedron *next = Q->next;
	    Q->next = NULL;

	    f = barvinok_enumerate_ev_f(Q, C, options);
	    emul(f, eres);
	    free_evalue_refs(f);
	    free(f);

	    Q->next = next;
	}
    }

    goto out;
}

evalue* barvinok_enumerate_with_options(Polyhedron *P, Polyhedron* C,
					struct barvinok_options *options)
{
    Polyhedron *next, *Cnext, *CA;
    Polyhedron *Porig = P;
    evalue *eres;

    if (P->next)
	fprintf(stderr,
    "barvinok_enumerate: input is a union; only first polyhedron is enumerated\n");

    if (C->next)
	fprintf(stderr,
    "barvinok_enumerate: context is a union; only first polyhedron is considered\n");

    Cnext = C->next;
    C->next = NULL;
    CA = align_context(C, P->Dimension, options->MaxRays);
    next = P->next;
    P->next = NULL;
    P = DomainIntersection(P, CA, options->MaxRays);
    Porig->next = next;
    Polyhedron_Free(CA);

    eres = enumerate(P, C, options);

    C->next = Cnext;

    return eres;
}

evalue* barvinok_enumerate_ev(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    evalue *E;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    E = barvinok_enumerate_with_options(P, C, options);
    barvinok_options_free(options);
    return E;
}

evalue *Param_Polyhedron_Enumerate(Param_Polyhedron *PP, Polyhedron *P,
				   Polyhedron *C,
				   struct barvinok_options *options)
{
    evalue *eres;
    Param_Domain *D;
    unsigned nparam = C->Dimension;
    unsigned dim = P->Dimension - nparam;

    ALLOC(evalue, eres);
    value_init(eres->d);
    value_set_si(eres->d, 0);

    int nd;
    for (nd = 0, D=PP->D; D; ++nd, D=D->next);
    struct section { Polyhedron *D; evalue E; };
    section *s = new section[nd];

    enumerator_base *et = NULL;
try_again:
    if (et)
	delete et;

    et = enumerator_base::create(P, dim, PP->nbV, options);

    Polyhedron *TC = true_context(P, C, options->MaxRays);
    FORALL_REDUCED_DOMAIN(PP, TC, nd, options, i, D, rVD)
	Param_Vertices *V;

	value_init(s[i].E.d);
	evalue_set_si(&s[i].E, 0, 1);
	s[i].D = rVD;

	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    if (!et->vE[_i])
		try {
		    et->decompose_at(V, _i, options);
		} catch (OrthogonalException &e) {
		    FORALL_REDUCED_DOMAIN_RESET;
		    for (; i >= 0; --i) {
			free_evalue_refs(&s[i].E);
			Domain_Free(s[i].D);
		    }
		    goto try_again;
		}
	    eadd(et->vE[_i] , &s[i].E);
	END_FORALL_PVertex_in_ParamPolyhedron;
	evalue_range_reduction_in_domain(&s[i].E, rVD);
    END_FORALL_REDUCED_DOMAIN
    Polyhedron_Free(TC);

    delete et;
    if (nd == 0)
	evalue_set_si(eres, 0, 1);
    else {
	eres->x.p = new_enode(partition, 2*nd, C->Dimension);
	for (int j = 0; j < nd; ++j) {
	    EVALUE_SET_DOMAIN(eres->x.p->arr[2*j], s[j].D);
	    value_clear(eres->x.p->arr[2*j+1].d);
	    eres->x.p->arr[2*j+1] = s[j].E;
	}
    }
    delete [] s;

    return eres;
}

static evalue* barvinok_enumerate_ev_f(Polyhedron *P, Polyhedron* C, 
				       barvinok_options *options)
{
    unsigned nparam = C->Dimension;
    bool do_scale = options->approximation_method == BV_APPROX_SCALE;

    if (options->approximation_method == BV_APPROX_VOLUME)
	return Param_Polyhedron_Volume(P, C, options);

    if (P->Dimension - nparam == 1 && !do_scale)
	return ParamLine_Length(P, C, options);

    Param_Polyhedron *PP = NULL;
    evalue *eres;

    if (do_scale) {
	eres = scale_bound(P, C, options);
	if (eres)
	    return eres;
    }

    PP = Polyhedron2Param_MR(P, C, options->MaxRays);

    if (do_scale)
	eres = scale(PP, P, C, options);
    else
	eres = Param_Polyhedron_Enumerate(PP, P, C, options);

    if (PP)
	Param_Polyhedron_Free(PP);

    return eres;
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

    Matrix *M = Matrix_Alloc(exist, exist);
    Vector_Copy(P->Constraint[r]+1+nvar, M->p[0], exist);
    Vector_Gcd(M->p[0], exist, &g);
    if (value_notone_p(g))
	Vector_AntiScale(M->p[0], M->p[0], g, exist);
    value_clear(g);

    int ok = unimodular_complete(M, 1);
    assert(ok);
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
		      unsigned exist, unsigned nparam, barvinok_options *options)
{
#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Or\n");
#endif /* DEBUG_ER */

    Polyhedron *N = D->next;
    D->next = 0;
    evalue *EP = 
	barvinok_enumerate_e_with_options(D, exist, nparam, options);
    Polyhedron_Free(D);

    for (D = N; D; D = N) {
	N = D->next;
	D->next = 0;

	evalue *EN = 
	    barvinok_enumerate_e_with_options(D, exist, nparam, options);

	eor(EN, EP);
	free_evalue_refs(EN); 
	free(EN);
	Polyhedron_Free(D);
    }

    reduce_evalue(EP);

    return EP;
}

static evalue* enumerate_sum(Polyhedron *P,
		      unsigned exist, unsigned nparam, barvinok_options *options)
{
    int nvar = P->Dimension - exist - nparam;
    int toswap = nvar < exist ? nvar : exist;
    for (int i = 0; i < toswap; ++i)
	SwapColumns(P, 1 + i, nvar+exist - i);
    nparam += nvar;

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Sum\n");
#endif /* DEBUG_ER */

    evalue *EP = barvinok_enumerate_e_with_options(P, exist, nparam, options);

    evalue_split_domains_into_orthants(EP, options->MaxRays);
    reduce_evalue(EP);
    evalue_range_reduction(EP);

    evalue_frac2floor2(EP, 1);

    evalue *sum = esum(EP, nvar);

    free_evalue_refs(EP); 
    free(EP);
    EP = sum;

    evalue_range_reduction(EP);

    return EP;
}

static evalue* split_sure(Polyhedron *P, Polyhedron *S,
		      unsigned exist, unsigned nparam, barvinok_options *options)
{
    int nvar = P->Dimension - exist - nparam;

    Matrix *M = Matrix_Alloc(exist, S->Dimension+2);
    for (int i = 0; i < exist; ++i)
	value_set_si(M->p[i][nvar+i+1], 1);
    Polyhedron *O = S;
    S = DomainAddRays(S, M, options->MaxRays);
    Polyhedron_Free(O);
    Polyhedron *F = DomainAddRays(P, M, options->MaxRays);
    Polyhedron *D = DomainDifference(F, S, options->MaxRays);
    O = D;
    D = Disjoint_Domain(D, 0, options->MaxRays);
    Polyhedron_Free(F);
    Domain_Free(O);
    Matrix_Free(M);

    M = Matrix_Alloc(P->Dimension+1-exist, P->Dimension+1);
    for (int j = 0; j < nvar; ++j)
	value_set_si(M->p[j][j], 1);
    for (int j = 0; j < nparam+1; ++j)
	value_set_si(M->p[nvar+j][nvar+exist+j], 1);
    Polyhedron *T = Polyhedron_Image(S, M, options->MaxRays);
    evalue *EP = barvinok_enumerate_e_with_options(T, 0, nparam, options);
    Polyhedron_Free(S);
    Polyhedron_Free(T);
    Matrix_Free(M);

    for (Polyhedron *Q = D; Q; Q = Q->next) {
	Polyhedron *N = Q->next;
	Q->next = 0;
	T = DomainIntersection(P, Q, options->MaxRays);
	evalue *E = barvinok_enumerate_e_with_options(T, exist, nparam, options);
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
		      unsigned exist, unsigned nparam, barvinok_options *options)
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
	S = AddConstraints(M->p[0], c, S, options->MaxRays);
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

    return split_sure(P, S, exist, nparam, options);
}

static evalue* enumerate_sure2(Polyhedron *P,
		      unsigned exist, unsigned nparam, barvinok_options *options)
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
    Polyhedron * F = Rays2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);

    Polyhedron *I = DomainIntersection(F, P, options->MaxRays);
    Polyhedron_Free(F);

#ifdef DEBUG_ER
    fprintf(stderr, "\nER: Sure2\n");
#endif /* DEBUG_ER */

    return split_sure(P, I, exist, nparam, options);
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
		      unsigned exist, unsigned nparam, barvinok_options *options)
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
    Polyhedron *S = AddConstraints(M->p[0], 2, P, options->MaxRays);
    evalue *EP = barvinok_enumerate_e_with_options(S, exist, nparam, options);
    Polyhedron_Free(S);
    Matrix_Free(M);

    return enumerate_cyclic(P, exist, nparam, EP, 0, i, options->MaxRays);
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
		      unsigned exist, unsigned nparam, barvinok_options *options)
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

    P = Rays2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    evalue *EP = barvinok_enumerate_e_with_options(P, exist, nparam, options);
    Polyhedron_Free(P);
    value_clear(one);

    return EP;
}

static evalue* enumerate_redundant_ray(Polyhedron *P,
		      unsigned exist, unsigned nparam, barvinok_options *options)
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
		return enumerate_remove_ray(P, r, exist, nparam, options);
	    }

	    value_division(m, P->Ray[r2][1+nvar+exist+i1], 
			      P->Ray[r][1+nvar+exist+i1]);
	    value_multiply(m, m, P->Ray[r][1+nvar+exist+i1]);
	    /* r divides r2 => r2 redundant */
	    if (value_eq(m, P->Ray[r2][1+nvar+exist+i1])) {
		value_clear(m);
		return enumerate_remove_ray(P, r2, exist, nparam, options);
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
		      unsigned exist, unsigned nparam, barvinok_options *options)
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
	    return enumerate_sum(P, exist, nparam, options);
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
    Polyhedron *S = Rays2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    Polyhedron *D = DomainDifference(P, S, options->MaxRays);
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
    D = AddConstraints(M->p[0], 2, P, options->MaxRays);
    // Polyhedron_Print(stderr, P_VALUE_FMT, D);
    value_subtract(M->p[0][1+P->Dimension], M->p[0][1+P->Dimension], 
		    P->Ray[r][1+nvar+exist+i]);
    // Matrix_Print(stderr, P_VALUE_FMT, M);
    S = AddConstraints(M->p[0], 1, P, options->MaxRays);
    // Polyhedron_Print(stderr, P_VALUE_FMT, S);
    Matrix_Free(M);

    evalue *EP = barvinok_enumerate_e_with_options(D, exist, nparam, options);
    Polyhedron_Free(D);
    value_clear(one);
    value_clear(m);

    if (value_notone_p(P->Ray[r][1+nvar+exist+i]))
	EP = enumerate_cyclic(P, exist, nparam, EP, r, i, options->MaxRays);
    else {
	M = Matrix_Alloc(1, nparam+2);
	value_set_si(M->p[0][0], 1);
	value_set_si(M->p[0][1+i], 1);
	enumerate_vd_add_ray(EP, M, options->MaxRays);
	Matrix_Free(M);
    }

    if (!emptyQ(S)) {
	evalue *E = barvinok_enumerate_e_with_options(S, exist, nparam, options);
	eadd(E, EP);
	free_evalue_refs(E);
	free(E);
    }
    Polyhedron_Free(S);

    if (R) {
	assert(nvar == 0);
	evalue *ER = enumerate_or(R, exist, nparam, options);
	eor(ER, EP);
	free_evalue_refs(ER);
	free(ER);
    }

    return EP;
}

static evalue* enumerate_vd(Polyhedron **PA,
		      unsigned exist, unsigned nparam, barvinok_options *options)
{
    Polyhedron *P = *PA;
    int nvar = P->Dimension - exist - nparam;
    Param_Polyhedron *PP = NULL;
    Polyhedron *C = Universe_Polyhedron(nparam);
    Polyhedron *CEq;
    Matrix *CT;
    Polyhedron *PR = P;
    PP = Polyhedron2Param_Domain(PR,C, options->MaxRays);
    Polyhedron_Free(C);

    int nd;
    Param_Domain *D, *last;
    Value c;
    value_init(c);
    for (nd = 0, D=PP->D; D; D=D->next, ++nd)
	;

    Polyhedron **VD = new Polyhedron_p[nd];
    Polyhedron *TC = true_context(P, C, options->MaxRays);
    FORALL_REDUCED_DOMAIN(PP, TC, nd, options, i, D, rVD)
	VD[nd++] = rVD;
	last = D;
    END_FORALL_REDUCED_DOMAIN
    Polyhedron_Free(TC);

    evalue *EP = 0;

    if (nd == 0)
	EP = evalue_zero();

    /* This doesn't seem to have any effect */
    if (nd == 1) {
	Polyhedron *CA = align_context(VD[0], P->Dimension, options->MaxRays);
	Polyhedron *O = P;
	P = DomainIntersection(P, CA, options->MaxRays);
	if (O != *PA)
	    Polyhedron_Free(O);
	Polyhedron_Free(CA);
	if (emptyQ(P))
	    EP = evalue_zero();
    }

    if (PR != *PA)
	Polyhedron_Free(PR);
    PR = 0;

    if (!EP && nd > 1) {
#ifdef DEBUG_ER
	fprintf(stderr, "\nER: VD\n");
#endif /* DEBUG_ER */
	for (int i = 0; i < nd; ++i) {
	    Polyhedron *CA = align_context(VD[i], P->Dimension, options->MaxRays);
	    Polyhedron *I = DomainIntersection(P, CA, options->MaxRays);

	    if (i == 0)
		EP = barvinok_enumerate_e_with_options(I, exist, nparam, options);
	    else {
		evalue *E = barvinok_enumerate_e_with_options(I, exist, nparam,
							      options);
		eadd(E, EP);
		free_evalue_refs(E); 
		free(E);
	    }
	    Polyhedron_Free(I);
	    Polyhedron_Free(CA);
	}
    }

    for (int i = 0; i < nd; ++i)
	Polyhedron_Free(VD[i]);
    delete [] VD;
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
						      options->MaxRays);
		    if (emptyQ(para)) {
			Polyhedron_Free(para);
			continue;
		    }
		    Polyhedron *pos, *neg;
		    value_set_si(M->p[0][0], 1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    neg = AddConstraints(M->p[0], 1, P, options->MaxRays);
		    value_set_si(f, -1);
		    Vector_Scale(M->p[0]+1, M->p[0]+1, f, 
				 P->Dimension+1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, options->MaxRays);
		    if (emptyQ(neg) && emptyQ(pos)) {
			Polyhedron_Free(para);
			Polyhedron_Free(pos);
			Polyhedron_Free(neg);
			continue;
		    }
#ifdef DEBUG_ER
		    fprintf(stderr, "\nER: Order\n");
#endif /* DEBUG_ER */
		    EP = barvinok_enumerate_e_with_options(para, exist, nparam,
							   options);
		    evalue *E;
		    if (!emptyQ(pos)) {
			E = barvinok_enumerate_e_with_options(pos, exist, nparam,
								options);
			eadd(E, EP);
			free_evalue_refs(E); 
			free(E);
		    }
		    if (!emptyQ(neg)) {
			E = barvinok_enumerate_e_with_options(neg, exist, nparam,
								options);
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
		    neg = AddConstraints(M->p[0], 1, P, options->MaxRays);
		    value_set_si(f, -1);
		    Vector_Scale(M->p[0]+1, M->p[0]+1, f, 
				 P->Dimension+1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, options->MaxRays);
		    if (emptyQ(neg) || emptyQ(pos)) {
			Polyhedron_Free(pos);
			Polyhedron_Free(neg);
			continue;
		    }
		    Polyhedron_Free(pos);
		    value_increment(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, options->MaxRays);
#ifdef DEBUG_ER
		    fprintf(stderr, "\nER: Vertex\n");
#endif /* DEBUG_ER */
		    pos->next = neg;
		    EP = enumerate_or(pos, exist, nparam, options);
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
		    neg = AddConstraints(M->p[0], 1, P, options->MaxRays);
		    value_set_si(f, -1);
		    Vector_Scale(M->p[0]+1, M->p[0]+1, f, 
				 P->Dimension+1);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    value_decrement(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, options->MaxRays);
		    if (emptyQ(neg) || emptyQ(pos)) {
			Polyhedron_Free(pos);
			Polyhedron_Free(neg);
			continue;
		    }
		    Polyhedron_Free(pos);
		    value_increment(M->p[0][P->Dimension+1],
				    M->p[0][P->Dimension+1]);
		    pos = AddConstraints(M->p[0], 1, P, options->MaxRays);
#ifdef DEBUG_ER
		    fprintf(stderr, "\nER: ParamVertex\n");
#endif /* DEBUG_ER */
		    pos->next = neg;
		    EP = enumerate_or(pos, exist, nparam, options);
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

evalue* barvinok_enumerate_pip(Polyhedron *P, unsigned exist, unsigned nparam,
			       unsigned MaxRays)
{
    evalue *E;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    E = barvinok_enumerate_pip_with_options(P, exist, nparam, options);
    barvinok_options_free(options);
    return E;
}

#ifndef HAVE_PIPLIB
evalue *barvinok_enumerate_pip_with_options(Polyhedron *P,
		  unsigned exist, unsigned nparam, struct barvinok_options *options)
{
    return 0;
}
#else
evalue *barvinok_enumerate_pip_with_options(Polyhedron *P,
		  unsigned exist, unsigned nparam, struct barvinok_options *options)
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
	E = barvinok_enumerate_e_with_options(Q, exist, nparam, options);
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
		      unsigned exist, unsigned nparam, barvinok_options *options);

#ifdef DEBUG_ER
static int er_level = 0;

evalue* barvinok_enumerate_e_with_options(Polyhedron *P, 
			  unsigned exist, unsigned nparam, barvinok_options *options)
{
    fprintf(stderr, "\nER: level %i\n", er_level);

    Polyhedron_PrintConstraints(stderr, P_VALUE_FMT, P);
    fprintf(stderr, "\nE %d\nP %d\n", exist, nparam);
    ++er_level;
    P = DomainConstraintSimplify(Polyhedron_Copy(P), options->MaxRays);
    evalue *EP = barvinok_enumerate_e_r(P, exist, nparam, options);
    Polyhedron_Free(P);
    --er_level;
    return EP;
}
#else
evalue* barvinok_enumerate_e_with_options(Polyhedron *P, 
			  unsigned exist, unsigned nparam, barvinok_options *options)
{
    P = DomainConstraintSimplify(Polyhedron_Copy(P), options->MaxRays);
    evalue *EP = barvinok_enumerate_e_r(P, exist, nparam, options);
    Polyhedron_Free(P);
    return EP;
}
#endif

evalue* barvinok_enumerate_e(Polyhedron *P, unsigned exist, unsigned nparam,
			     unsigned MaxRays)
{
    evalue *E;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    E = barvinok_enumerate_e_with_options(P, exist, nparam, options);
    barvinok_options_free(options);
    return E;
}

static evalue* barvinok_enumerate_e_r(Polyhedron *P, 
			  unsigned exist, unsigned nparam, barvinok_options *options)
{
    if (exist == 0) {
	Polyhedron *U = Universe_Polyhedron(nparam);
	evalue *EP = barvinok_enumerate_with_options(P, U, options);
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
	barvinok_count_with_options(P, &EP->x.n, options);
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
	    Polyhedron *T = rotate_along(P, r, nvar, exist, options->MaxRays);
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Equality\n");
#endif /* DEBUG_ER */
	    evalue *EP = barvinok_enumerate_e_with_options(T, exist-1, nparam,
							   options);
	    Polyhedron_Free(T);
	    return EP;
	} else {
#ifdef DEBUG_ER
	    fprintf(stderr, "\nER: Fixed\n");
#endif /* DEBUG_ER */
	    if (first == 0)
		return barvinok_enumerate_e_with_options(P, exist-1, nparam,
							 options);
	    else {
		Polyhedron *T = Polyhedron_Copy(P);
		SwapColumns(T, nvar+1, nvar+1+first);
		evalue *EP = barvinok_enumerate_e_with_options(T, exist-1, nparam,
							       options);
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
		    Polyhedron *T = AddConstraints(row->p, 1, P, options->MaxRays);
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
			Polyhedron *T = AddConstraints(row->p, 1, P,
						       options->MaxRays);
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
			    Polyhedron *T = AddConstraints(row->p, 1, P,
							   options->MaxRays);
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
	    Polyhedron *T = Polyhedron_Image(P, M, options->MaxRays);
	    Matrix_Free(M);
	    evalue *EP = barvinok_enumerate_e_with_options(T, exist-1, nparam,
							   options);
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
		return barvinok_enumerate_e_with_options(P, exist-1, nparam,
							 options);
	    else {
		Polyhedron *T = Polyhedron_Copy(P);
		SwapColumns(T, nvar+1, nvar+1+i);
		evalue *EP = barvinok_enumerate_e_with_options(T, exist-1, nparam,
							       options);
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
	    Polyhedron *T = rotate_along(P, r, nvar, exist, options->MaxRays);
	    evalue *EP = barvinok_enumerate_e_with_options(T, exist-1, nparam,
							   options);
	    Polyhedron_Free(T);
	    return EP;
	}
    for (int i = 0; i < exist; ++i)
	if (info[i] & INDEPENDENT) {
	    Polyhedron *pos, *neg;

	    /* Find constraint again and split off negative part */

	    if (SplitOnVar(P, i, nvar, exist, options->MaxRays,
			   row, f, true, &pos, &neg)) {
#ifdef DEBUG_ER
		fprintf(stderr, "\nER: Split\n");
#endif /* DEBUG_ER */

		evalue *EP = 
		    barvinok_enumerate_e_with_options(neg, exist-1, nparam, options);
		evalue *E = 
		    barvinok_enumerate_e_with_options(pos, exist, nparam, options);
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

    EP = enumerate_line(P, exist, nparam, options);
    if (EP)
	goto out;

    EP = barvinok_enumerate_pip_with_options(P, exist, nparam, options);
    if (EP)
	goto out;

    EP = enumerate_redundant_ray(P, exist, nparam, options);
    if (EP)
	goto out;

    EP = enumerate_sure(P, exist, nparam, options);
    if (EP)
	goto out;

    EP = enumerate_ray(P, exist, nparam, options);
    if (EP)
	goto out;

    EP = enumerate_sure2(P, exist, nparam, options);
    if (EP)
	goto out;

    F = unfringe(P, options->MaxRays);
    if (!PolyhedronIncludes(F, P)) {
#ifdef DEBUG_ER
	fprintf(stderr, "\nER: Fringed\n");
#endif /* DEBUG_ER */
	EP = barvinok_enumerate_e_with_options(F, exist, nparam, options);
	Polyhedron_Free(F);
	goto out;
    }
    Polyhedron_Free(F);

    if (nparam)
	EP = enumerate_vd(&P, exist, nparam, options);
    if (EP)
	goto out2;

    if (nvar != 0) {
	EP = enumerate_sum(P, exist, nparam, options);
	goto out2;
    }

    assert(nvar == 0);

    int i;
    Polyhedron *pos, *neg;
    for (i = 0; i < exist; ++i)
	if (SplitOnVar(P, i, nvar, exist, options->MaxRays,
		       row, f, false, &pos, &neg))
	    break;

    assert (i < exist);

    pos->next = neg;
    EP = enumerate_or(pos, exist, nparam, options);

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
static Polyhedron *remove_more_equalities(Polyhedron *P, unsigned nparam,
					  Matrix **CP, unsigned MaxRays)
{
    Polyhedron *Q = P;
    remove_all_equalities(&P, NULL, CP, NULL, nparam, MaxRays);
    if (P != Q)
	Polyhedron_Free(Q);
    return P;
}

/* frees P */
static gen_fun *series(Polyhedron *P, unsigned nparam, barvinok_options *options)
{
    Matrix *CP = NULL;
    gen_fun *gf;

    if (emptyQ2(P)) {
	Polyhedron_Free(P);
	return new gen_fun;
    }

    assert(!Polyhedron_is_unbounded(P, nparam, options->MaxRays));
    assert(P->NbBid == 0);
    assert(Polyhedron_has_revlex_positive_rays(P, nparam));
    if (P->NbEq != 0)
	P = remove_more_equalities(P, nparam, &CP, options->MaxRays);
    assert(P->NbEq == 0);
    if (CP)
	nparam = CP->NbColumns-1;

    if (nparam == 0) {
	Value c;
	value_init(c);
	barvinok_count_with_options(P, &c, options);
	gf = new gen_fun(c);
	value_clear(c);
    } else {
	gf_base *red;
	red = gf_base::create(Polyhedron_Project(P, nparam),
			      P->Dimension, nparam, options);
	POL_ENSURE_VERTICES(P);
	red->start_gf(P, options);
	gf = red->gf;
	delete red;
    }
    if (CP) {
	gf->substitute(CP);
	Matrix_Free(CP);
    }
    Polyhedron_Free(P);
    return gf;
}

gen_fun * barvinok_series_with_options(Polyhedron *P, Polyhedron* C,
				       barvinok_options *options)
{
    Polyhedron *CA;
    unsigned nparam = C->Dimension;
    gen_fun *gf;

    CA = align_context(C, P->Dimension, options->MaxRays);
    P = DomainIntersection(P, CA, options->MaxRays);
    Polyhedron_Free(CA);

    gf = series(P, nparam, options);

    return gf;
}

gen_fun * barvinok_series(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    gen_fun *gf;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    gf = barvinok_series_with_options(P, C, options);
    barvinok_options_free(options);
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
	assert(!Polyhedron_is_unbounded(P, nparam, MaxRays));
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

gen_fun* barvinok_enumerate_union_series_with_options(Polyhedron *D, Polyhedron* C, 
						      barvinok_options *options)
{
    Polyhedron *conv, *D2;
    Polyhedron *CA;
    gen_fun *gf = NULL, *gf2;
    unsigned nparam = C->Dimension;
    ZZ one, mone;
    one = 1;
    mone = -1;

    CA = align_context(C, D->Dimension, options->MaxRays);
    D = DomainIntersection(D, CA, options->MaxRays);
    Polyhedron_Free(CA);

    D2 = skew_into_positive_orthant(D, nparam, options->MaxRays);
    for (Polyhedron *P = D2; P; P = P->next) {
	assert(P->Dimension == D2->Dimension);
	gen_fun *P_gf;

	P_gf = series(Polyhedron_Copy(P), nparam, options);
	if (!gf)
	    gf = P_gf;
	else {
	    gf->add_union(P_gf, options);
	    delete P_gf;
	}
    }
    /* we actually only need the convex union of the parameter space
     * but the reducer classes currently expect a polyhedron in
     * the combined space
     */
    Polyhedron_Free(gf->context);
    gf->context = DomainConvex(D2, options->MaxRays);

    gf2 = gf->summate(D2->Dimension - nparam, options);

    delete gf;
    if (D != D2)
	Domain_Free(D2);
    Domain_Free(D);
    return gf2;
}

gen_fun* barvinok_enumerate_union_series(Polyhedron *D, Polyhedron* C, 
					 unsigned MaxRays)
{
    gen_fun *gf;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    gf = barvinok_enumerate_union_series_with_options(D, C, options);
    barvinok_options_free(options);
    return gf;
}

evalue* barvinok_enumerate_union(Polyhedron *D, Polyhedron* C, unsigned MaxRays)
{
    evalue *EP;
    gen_fun *gf = barvinok_enumerate_union_series(D, C, MaxRays);
    EP = *gf;
    delete gf;
    return EP;
}
