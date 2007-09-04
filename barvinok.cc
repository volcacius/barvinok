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
#include "config.h"
#include <barvinok/barvinok.h>
#include <barvinok/genfun.h>
#include <barvinok/options.h>
#include <barvinok/sample.h>
#include "bfcounter.h"
#include "conversion.h"
#include "counter.h"
#include "tcounter.h"
#include "decomposer.h"
#include "lattice_point.h"
#include "reduce_domain.h"
#include "remove_equalities.h"
#include "scale.h"
#include "volume.h"
#include "bernoulli.h"
#include "param_util.h"

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
    dpoly_n(int d) {
	Value d0, one;
	value_init(d0);
	value_init(one);
	value_set_si(one, 1);
	coeff = Matrix_Alloc(d+1, d+1+1);
	value_set_si(coeff->p[0][0], 1);
	value_set_si(coeff->p[0][d+1], 1);
	for (int i = 1; i <= d; ++i) {
	    value_multiply(coeff->p[i][0], coeff->p[i-1][0], d0);
	    Vector_Combine(coeff->p[i-1], coeff->p[i-1]+1, coeff->p[i]+1,
			   one, d0, i);
	    value_set_si(coeff->p[i][d+1], i);
	    value_multiply(coeff->p[i][d+1], coeff->p[i][d+1], coeff->p[i-1][d+1]);
	    value_decrement(d0, d0);
	}
	value_clear(d0);
	value_clear(one);
    }
    void div(dpoly& d, Vector *count, int sign) {
	int len = coeff->NbRows;
	Matrix * c = Matrix_Alloc(coeff->NbRows, coeff->NbColumns);
	Value tmp;
	value_init(tmp);
	for (int i = 0; i < len; ++i) {
	    Vector_Copy(coeff->p[i], c->p[i], len+1);
	    for (int j = 1; j <= i; ++j) {
		value_multiply(tmp, d.coeff->p[j], c->p[i][len]);
		value_oppose(tmp, tmp);
		Vector_Combine(c->p[i], c->p[i-j], c->p[i],
			       c->p[i-j][len], tmp, len);
		value_multiply(c->p[i][len], c->p[i][len], c->p[i-j][len]);
	    }
	    value_multiply(c->p[i][len], c->p[i][len], d.coeff->p[0]);
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

static void print_int_vector(int *v, int len, const char *name)
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
           const char * test[] = {"a", "b"};
           print_evalue(stderr, bfet->factors[j], test);
           fprintf(stderr, "\n");
	}
    }
}

struct bfcounter : public bfcounter_base {
    mpq_t count;
    Value tz;

    bfcounter(unsigned dim) : bfcounter_base(dim) {
	mpq_init(count);
	lower = 1;
	value_init(tz);
    }
    ~bfcounter() {
	mpq_clear(count);
	value_clear(tz);
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

	zz2value(factors[j][0], tz);
	dpoly D(total_power, tz, 1);
	for (int k = 1; k < v[i]->powers[j]; ++k) {
	    zz2value(factors[j][0], tz);
	    dpoly fact(total_power, tz, 1);
	    D *= fact;
	}
	for ( ; ++j < nf; )
	    for (int k = 0; k < v[i]->powers[j]; ++k) {
		zz2value(factors[j][0], tz);
		dpoly fact(total_power, tz, 1);
		D *= fact;
	    }

	for (int k = 0; k < v[i]->terms.NumRows(); ++k) {
	    zz2value(v[i]->terms[k][0], tz);
	    dpoly n(total_power, tz);
	    mpq_set_si(tcount, 0, 1);
	    n.div(D, tcount, 1);
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
    if (options->incremental_specialization == BV_SPECIALIZATION_BF)
	cnt = new bfcounter(P->Dimension);
    else if (options->incremental_specialization == BV_SPECIALIZATION_DF)
	cnt = new icounter(P->Dimension);
    else if (options->incremental_specialization == BV_SPECIALIZATION_TODD)
	cnt = new tcounter(P->Dimension, options->max_index);
    else
	cnt = new counter(P->Dimension, options->max_index);
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
    term_info num;
    Vector *c;
    mpq_t count;
    Value tz;

    enumerator(Polyhedron *P, unsigned dim, unsigned nbV) :
		vertex_decomposer(P, nbV, *this), enumerator_base(dim, this) {
	this->P = P;
	this->nbV = nbV;
	randomvector(P, lambda, dim);
	den.SetLength(dim);
	c = Vector_Alloc(dim+2);

	mpq_init(count);
	value_init(tz);
    }

    ~enumerator() {
	mpq_clear(count);
	Vector_Free(c);
	value_clear(tz);
    }

    virtual void handle(const signed_cone& sc, barvinok_options *options);
};

void enumerator::handle(const signed_cone& sc, barvinok_options *options)
{
    int sign = sc.sign;
    int r = 0;
    assert(sc.rays.NumRows() == dim);
    for (int k = 0; k < dim; ++k) {
	if (lambda * sc.rays[k] == 0)
	    throw Orthogonal;
    }

    lattice_point(V, sc.rays, lambda, &num, sc.det, sc.closed, options);
    den = sc.rays * lambda;

    if (dim % 2)
	sign = -sign;

    zz2value(den[0], tz);
    dpoly n(dim, tz, 1);
    for (int k = 1; k < dim; ++k) {
	zz2value(den[k], tz);
	dpoly fact(dim, tz, 1);
	n *= fact;
    }
    if (num.E != NULL) {
	dpoly_n d(dim);
	d.div(n, c, sign);
	for (unsigned long i = 0; i < sc.det; ++i) {
	    evalue *EV = evalue_polynomial(c, num.E[i]);
	    eadd(EV, vE[vert]);
	    free_evalue_refs(EV);
	    free(EV);
	    free_evalue_refs(num.E[i]);
	    delete num.E[i];
	}
	delete [] num.E; 
    } else {
	mpq_set_si(count, 0, 1);
	if (num.constant.length() == 1) {
	    zz2value(num.constant[0], tz);
	    dpoly d(dim, tz);
	    d.div(n, count, sign);
	} else {
	    dpoly_n d(dim);
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
    Value tz;

    ienumerator(Polyhedron *P, unsigned dim, unsigned nbV) :
		vertex_decomposer(P, nbV, *this), ienumerator_base(dim, this) {
	vertex.SetDims(1, dim);

	den.SetDims(dim, dim);
	mpq_init(tcount);
	value_init(tz);
    }

    ~ienumerator() {
	mpq_clear(tcount);
	value_clear(tz);
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

	zz2value(num_s[0], tz);
	dpoly n(no_param, tz);
	zz2value(den_s[k], tz);
	dpoly D(no_param, tz, 1);
	for ( ; ++k < len; )
	    if (den_p[k] == 0) {
		zz2value(den_s[k], tz);
		dpoly fact(no_param, tz, 1);
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
		n.div(D, tcount, 1);

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

		zz2value(den_s[k], tz);
		dpoly pd(no_param-1, tz, 1);

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

static evalue* barvinok_enumerate_ev_f(Polyhedron *P, Polyhedron* C, 
				       barvinok_options *options);

/* Destroys C */
static evalue* barvinok_enumerate_cst(Polyhedron *P, Polyhedron* C, 
				      struct barvinok_options *options)
{
    evalue *eres;

    if (emptyQ2(C)) {
	Polyhedron_Free(C);
	return evalue_zero();
    }

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

static evalue* enumerate(Polyhedron *P, Polyhedron* C,
					struct barvinok_options *options)
{
    Polyhedron *next;
    Polyhedron *Porig = P;
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

    if (C->Dimension == 0 || emptyQ(P) || emptyQ(C)) {
constant:
	if (CEq == Porig)
	    CEq = Polyhedron_Copy(CEq);
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
	if (P != Porig)
	    Domain_Free(P);
	if (C != Corig)
	    Polyhedron_Free(C);
	   
	return eres;
    }
    if (Polyhedron_is_unbounded(P, nparam, options->MaxRays))
	goto constant;

    if (P->NbEq != 0) {
	Matrix *f;
	P = remove_equalities_p(Polyhedron_Copy(P), P->Dimension-nparam, &f,
				options->MaxRays);
	mask(f, &factor, options);
	Matrix_Free(f);
    }
    if (P->Dimension == nparam) {
	CEq = P;
	P = Universe_Polyhedron(0);
	goto constant;
    }
    if (P->NbEq != 0 || C->NbEq != 0) {
	Polyhedron *Q = P;
	Polyhedron *D = C;
	remove_all_equalities(&P, &C, &CP, NULL, nparam, options->MaxRays);
	if (C != D && D != Corig)
	    Polyhedron_Free(D);
	if (P != Q && Q != Porig)
	    Domain_Free(Q);
	eres = enumerate(P, C, options);
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
	if (P != Porig)
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
    Polyhedron *next, *Cnext, *C1;
    Polyhedron *Corig = C;
    evalue *eres;

    if (P->next)
	fprintf(stderr,
    "barvinok_enumerate: input is a union; only first polyhedron is enumerated\n");

    if (C->next)
	fprintf(stderr,
    "barvinok_enumerate: context is a union; only first polyhedron is considered\n");

    Cnext = C->next;
    C->next = NULL;
    C1 = Polyhedron_Project(P, C->Dimension);
    C = DomainIntersection(C, C1, options->MaxRays);
    Polyhedron_Free(C1);
    next = P->next;
    P->next = NULL;

    if (options->approximation_method == BV_APPROX_BERNOULLI)
	eres = Bernoulli_sum(P, C, options);
    else
	eres = enumerate(P, C, options);
    Domain_Free(C);

    P->next= next;
    Corig->next = Cnext;

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

    int nd;
    for (nd = 0, D=PP->D; D; ++nd, D=D->next);
    evalue_section *s = new evalue_section[nd];

    enumerator_base *et = NULL;
try_again:
    if (et)
	delete et;

    et = enumerator_base::create(P, dim, PP->nbV, options);

    Polyhedron *TC = true_context(P, C, options->MaxRays);
    FORALL_REDUCED_DOMAIN(PP, TC, nd, options, i, D, rVD)
	Param_Vertices *V;

	s[i].E = evalue_zero();
	s[i].D = rVD;

	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    if (!et->vE[_i])
		try {
		    et->decompose_at(V, _i, options);
		} catch (OrthogonalException &e) {
		    FORALL_REDUCED_DOMAIN_RESET;
		    for (; i >= 0; --i) {
			free_evalue_refs(s[i].E);
			free(s[i].E);
			Domain_Free(s[i].D);
		    }
		    goto try_again;
		}
	    eadd(et->vE[_i] , s[i].E);
	END_FORALL_PVertex_in_ParamPolyhedron;
	evalue_range_reduction_in_domain(s[i].E, rVD);
    END_FORALL_REDUCED_DOMAIN
    Polyhedron_Free(TC);

    delete et;
    eres = evalue_from_section_array(s, nd);
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

    PP = Polyhedron2Param_Polyhedron(P, C, options);

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

evalue* barvinok_enumerate_union(Polyhedron *D, Polyhedron* C, unsigned MaxRays)
{
    evalue *EP;
    gen_fun *gf = barvinok_enumerate_union_series(D, C, MaxRays);
    EP = *gf;
    delete gf;
    return EP;
}
