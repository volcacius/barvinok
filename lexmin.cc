#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <gmp.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/barvinok.h>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "decomposer.h"
#include "lattice_point.h"
#include "reduce_domain.h"
#include "mat_util.h"
#include "combine.h"
#include "edomain.h"
#include "evalue_util.h"
#include "sample.h"
#include "config.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

using std::vector;
using std::map;
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

/* RANGE : normal range for evalutations (-RANGE -> RANGE) */
#define RANGE 50

/* SRANGE : small range for evalutations */
#define SRANGE 15

/* if dimension >= BIDDIM, use SRANGE */
#define BIGDIM 5

/* VSRANGE : very small range for evalutations */
#define VSRANGE 5

/* if dimension >= VBIDDIM, use VSRANGE */
#define VBIGDIM 8

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
#define NO_EMPTINESS_CHECK  256
struct option lexmin_options[] = {
    { "verify",     no_argument,  0,  'T' },
    { "print-all",  no_argument,  0,  'A' },
    { "no-emptiness-check", no_argument, 0, NO_EMPTINESS_CHECK },
    { "min",   	    required_argument,  0,  'm' },
    { "max",   	    required_argument,  0,  'M' },
    { "range",      required_argument,  0,  'r' },
    { "version",    no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 : 
	  p->type == flooring ? 1 : 0;
}

struct indicator_term {
    int sign;
    int pos;
    mat_ZZ den;
    evalue **vertex;

    indicator_term(unsigned dim, int pos) {
	den.SetDims(0, dim);
	vertex = new evalue* [dim];
	this->pos = pos;
	sign = 0;
    }
    indicator_term(unsigned dim) {
	den.SetDims(dim, dim);
	vertex = new evalue* [dim];
	pos = -1;
    }
    indicator_term(const indicator_term& src) {
	sign = src.sign;
	pos = src.pos;
	den = src.den;
	unsigned dim = den.NumCols();
	vertex = new evalue* [dim];
	for (int i = 0; i < dim; ++i) {
	    vertex[i] = new evalue();
	    value_init(vertex[i]->d);
	    evalue_copy(vertex[i], src.vertex[i]);
	}
    }
    ~indicator_term() {
	unsigned dim = den.NumCols();
	for (int i = 0; i < dim; ++i) {
	    free_evalue_refs(vertex[i]);
	    delete vertex[i];
	}
	delete [] vertex;
    }
    void print(ostream& os, char **p) const;
    void substitute(Matrix *T);
    void normalize();
    void substitute(evalue *fract, evalue *val);
    void substitute(int pos, evalue *val);
    void reduce_in_domain(Polyhedron *D);
    bool is_opposite(indicator_term *neg);
};

bool indicator_term::is_opposite(indicator_term *neg)
{
    if (sign + neg->sign != 0)
	return false;
    if (den != neg->den)
	return false;
    for (int k = 0; k < den.NumCols(); ++k)
	if (!eequal(vertex[k], neg->vertex[k]))
	    return false;
    return true;
}

void indicator_term::reduce_in_domain(Polyhedron *D)
{
    for (int k = 0; k < den.NumCols(); ++k) {
	reduce_evalue_in_domain(vertex[k], D);
	if (evalue_range_reduction_in_domain(vertex[k], D))
	    reduce_evalue(vertex[k]);
    }
}

void indicator_term::print(ostream& os, char **p) const
{
    unsigned dim = den.NumCols();
    unsigned factors = den.NumRows();
    if (sign == 0)
	os << " s ";
    else if (sign > 0)
	os << " + ";
    else
	os << " - ";
    os << '[';
    for (int i = 0; i < dim; ++i) {
	if (i)
	    os << ',';
	evalue_print(os, vertex[i], p);
    }
    os << ']';
    for (int i = 0; i < factors; ++i) {
	os << " + t" << i << "*[";
	for (int j = 0; j < dim; ++j) {
	    if (j)
		os << ',';
	    os << den[i][j];
	}
	os << "]";
    }
    if (sign == 0)
	os << " (" << pos << ")";
}

/* Perform the substitution specified by T on the variables.
 * T has dimension (newdim+nparam+1) x (olddim + nparam + 1).
 * The substitution is performed as in gen_fun::substitute
 */
void indicator_term::substitute(Matrix *T)
{
    unsigned dim = den.NumCols();
    unsigned nparam = T->NbColumns - dim - 1;
    unsigned newdim = T->NbRows - nparam - 1;
    evalue **newvertex;
    mat_ZZ trans;
    matrix2zz(T, trans, newdim, dim);
    trans = transpose(trans);
    den *= trans;
    newvertex = new evalue* [newdim];

    vec_ZZ v;
    v.SetLength(nparam+1);

    evalue factor;
    value_init(factor.d);
    value_set_si(factor.d, 1);
    value_init(factor.x.n);
    for (int i = 0; i < newdim; ++i) {
	values2zz(T->p[i]+dim, v, nparam+1);
	newvertex[i] = multi_monom(v);

	for (int j = 0; j < dim; ++j) {
	    if (value_zero_p(T->p[i][j]))
		continue;
	    evalue term;
	    value_init(term.d);
	    evalue_copy(&term, vertex[j]);
	    value_assign(factor.x.n, T->p[i][j]);
	    emul(&factor, &term);
	    eadd(&term, newvertex[i]);
	    free_evalue_refs(&term);
	}
    }
    free_evalue_refs(&factor);
    for (int i = 0; i < dim; ++i) {
	free_evalue_refs(vertex[i]);
	delete vertex[i];
    }
    delete [] vertex;
    vertex = newvertex;
}

static void evalue_add_constant(evalue *e, ZZ v)
{
    Value tmp;
    value_init(tmp);

    /* go down to constant term */
    while (value_zero_p(e->d))
	e = &e->x.p->arr[type_offset(e->x.p)];
    /* and add v */
    zz2value(v, tmp);
    value_multiply(tmp, tmp, e->d);
    value_addto(e->x.n, e->x.n, tmp);

    value_clear(tmp);
}

/* Make all powers in denominator lexico-positive */
void indicator_term::normalize()
{
    vec_ZZ extra_vertex;
    extra_vertex.SetLength(den.NumCols());
    for (int r = 0; r < den.NumRows(); ++r) {
	for (int k = 0; k < den.NumCols(); ++k) {
	    if (den[r][k] == 0)
		continue;
	    if (den[r][k] > 0)
		break;
	    sign = -sign;
	    den[r] = -den[r];
	    extra_vertex += den[r];
	    break;
	}
    }
    for (int k = 0; k < extra_vertex.length(); ++k)
	if (extra_vertex[k] != 0)
	    evalue_add_constant(vertex[k], extra_vertex[k]);
}

static void substitute(evalue *e, evalue *fract, evalue *val)
{
    evalue *t;

    for (t = e; value_zero_p(t->d); t = &t->x.p->arr[type_offset(t->x.p)]) {
	if (t->x.p->type == fractional && eequal(&t->x.p->arr[0], fract))
	    break;
    }
    if (value_notzero_p(t->d))
	return;

    free_evalue_refs(&t->x.p->arr[0]);
    evalue *term = &t->x.p->arr[2];
    enode *p = t->x.p;
    value_clear(t->d);
    *t = t->x.p->arr[1];

    emul(val, term);
    eadd(term, e);
    free_evalue_refs(term);
    free(p);

    reduce_evalue(e);
}

void indicator_term::substitute(evalue *fract, evalue *val)
{
    unsigned dim = den.NumCols();
    for (int i = 0; i < dim; ++i) {
	::substitute(vertex[i], fract, val);
    }
}

static void substitute(evalue *e, int pos, evalue *val)
{
    evalue *t;

    for (t = e; value_zero_p(t->d); t = &t->x.p->arr[type_offset(t->x.p)]) {
	if (t->x.p->type == polynomial && t->x.p->pos == pos)
	    break;
    }
    if (value_notzero_p(t->d))
	return;

    evalue *term = &t->x.p->arr[1];
    enode *p = t->x.p;
    value_clear(t->d);
    *t = t->x.p->arr[0];

    emul(val, term);
    eadd(term, e);
    free_evalue_refs(term);
    free(p);

    reduce_evalue(e);
}

void indicator_term::substitute(int pos, evalue *val)
{
    unsigned dim = den.NumCols();
    for (int i = 0; i < dim; ++i) {
	::substitute(vertex[i], pos, val);
    }
}

struct indicator_constructor : public polar_decomposer, public vertex_decomposer {
    vec_ZZ vertex;
    vector<indicator_term*> *terms;
    Matrix *T;	/* Transformation to original space */
    Param_Polyhedron *PP;

    indicator_constructor(Polyhedron *P, unsigned dim, Param_Polyhedron *PP,
			  Matrix *T) :
		vertex_decomposer(P, PP->nbV, this), T(T), PP(PP) {
	vertex.SetLength(dim);
	terms = new vector<indicator_term*>[nbV];
    }
    ~indicator_constructor() {
	for (int i = 0; i < nbV; ++i)
	    for (int j = 0; j < terms[i].size(); ++j)
		delete terms[i][j];
	delete [] terms;
    }
    void substitute(Matrix *T);
    void normalize();
    void print(ostream& os, char **p);

    virtual void handle_polar(Polyhedron *P, int sign);
};

void indicator_constructor::handle_polar(Polyhedron *C, int s)
{
    unsigned dim = vertex.length();

    assert(C->NbRays-1 == dim);

    indicator_term *term = new indicator_term(dim);
    term->sign = s;
    terms[vert].push_back(term);

    lattice_point(V, C, vertex, term->vertex);

    for (int r = 0; r < dim; ++r) {
	values2zz(C->Ray[r]+1, term->den[r], dim);
	for (int j = 0; j < dim; ++j) {
	    if (term->den[r][j] == 0)
		continue;
	    if (term->den[r][j] > 0)
		break;
	    term->sign = -term->sign;
	    term->den[r] = -term->den[r];
	    vertex += term->den[r];
	    break;
	}
    }

    for (int i = 0; i < dim; ++i) {
	if (!term->vertex[i]) {
	    term->vertex[i] = new evalue();
	    value_init(term->vertex[i]->d);
	    value_init(term->vertex[i]->x.n);
	    zz2value(vertex[i], term->vertex[i]->x.n);
	    value_set_si(term->vertex[i]->d, 1);
	    continue;
	}
	if (vertex[i] == 0)
	    continue;
	evalue_add_constant(term->vertex[i], vertex[i]);
    }

    if (T) {
	term->substitute(T);
	term->normalize();
    }

    lex_order_rows(term->den);
}

void indicator_constructor::print(ostream& os, char **p)
{
    for (int i = 0; i < nbV; ++i)
	for (int j = 0; j < terms[i].size(); ++j) {
	    os << "i: " << i << ", j: " << j << endl;
	    terms[i][j]->print(os, p);
	    os << endl;
	}
}

void indicator_constructor::normalize()
{
    for (int i = 0; i < nbV; ++i)
	for (int j = 0; j < terms[i].size(); ++j) {
	    vec_ZZ vertex;
	    vertex.SetLength(terms[i][j]->den.NumCols());
	    for (int r = 0; r < terms[i][j]->den.NumRows(); ++r) {
		for (int k = 0; k < terms[i][j]->den.NumCols(); ++k) {
		    if (terms[i][j]->den[r][k] == 0)
			continue;
		    if (terms[i][j]->den[r][k] > 0)
			break;
		    terms[i][j]->sign = -terms[i][j]->sign;
		    terms[i][j]->den[r] = -terms[i][j]->den[r];
		    vertex += terms[i][j]->den[r];
		    break;
		}
	    }
	    lex_order_rows(terms[i][j]->den);
	    for (int k = 0; k < vertex.length(); ++k)
		if (vertex[k] != 0)
		    evalue_add_constant(terms[i][j]->vertex[k], vertex[k]);
	}
}

struct indicator;

enum order_sign { order_lt, order_le, order_eq, order_ge, order_gt, order_unknown };

struct partial_order {
    indicator *ind;

    map<indicator_term *, int > pred;
    map<indicator_term *, vector<indicator_term * > > lt;
    map<indicator_term *, vector<indicator_term * > > le;
    map<indicator_term *, vector<indicator_term * > > eq;

    map<indicator_term *, vector<indicator_term * > > pending;

    partial_order(indicator *ind) : ind(ind) {}
    void copy(const partial_order& order, 
	      map< indicator_term *, indicator_term * > old2new);

    order_sign compare(indicator_term *a, indicator_term *b);
    void set_equal(indicator_term *a, indicator_term *b);
    void unset_le(indicator_term *a, indicator_term *b);

    bool compared(indicator_term* a, indicator_term* b);
    void add(indicator_term* it, std::set<indicator_term *> *filter);
    void remove(indicator_term* it);

    void print(ostream& os, char **p);
};

void partial_order::unset_le(indicator_term *a, indicator_term *b)
{
    vector<indicator_term *>::iterator i;
    i = find(le[a].begin(), le[a].end(), b);
    le[a].erase(i);
    pred[b]--;
    i = find(pending[a].begin(), pending[a].end(), b);
    if (i != pending[a].end())
	pending[a].erase(i);
}

void partial_order::set_equal(indicator_term *a, indicator_term *b)
{
    if (eq[a].size() == 0)
	eq[a].push_back(a);
    if (eq[b].size() == 0)
	eq[b].push_back(b);
    a = eq[a][0];
    b = eq[b][0];
    assert(a != b);
    if (b < a) {
	indicator_term *c = a;
	a = b;
	b = c;
    }

    indicator_term *base = a;

    map<indicator_term *, vector<indicator_term * > >::iterator i;

    for (int j = 0; j < eq[b].size(); ++j) {
	eq[base].push_back(eq[b][j]);
	eq[eq[b][j]][0] = base;
    }
    eq[b].resize(1);

    i = lt.find(b);
    if (i != lt.end()) {
	for (int j = 0; j < lt[b].size(); ++j) {
	    if (find(eq[base].begin(), eq[base].end(), lt[b][j]) != eq[base].end())
		pred[lt[b][j]]--;
	    else if (find(lt[base].begin(), lt[base].end(), lt[b][j])
			!= lt[base].end())
		pred[lt[b][j]]--;
	    else
		lt[base].push_back(lt[b][j]);
	}
	lt.erase(i);
    }

    i = le.find(b);
    if (i != le.end()) {
	for (int j = 0; j < le[b].size(); ++j) {
	    if (find(eq[base].begin(), eq[base].end(), le[b][j]) != eq[base].end())
		pred[le[b][j]]--;
	    else if (find(le[base].begin(), le[base].end(), le[b][j])
			!= le[base].end())
		pred[le[b][j]]--;
	    else
		le[base].push_back(le[b][j]);
	}
	le.erase(i);
    }

    i = pending.find(base);
    if (i != pending.end()) {
	vector<indicator_term * > old = pending[base];
	pending[base].clear();
	for (int j = 0; j < old.size(); ++j) {
	    if (find(eq[base].begin(), eq[base].end(), old[j]) == eq[base].end())
		pending[base].push_back(old[j]);
	}
    }

    i = pending.find(b);
    if (i != pending.end()) {
	for (int j = 0; j < pending[b].size(); ++j) {
	    if (find(eq[base].begin(), eq[base].end(), pending[b][j]) == eq[base].end())
		pending[base].push_back(pending[b][j]);
	}
	pending.erase(i);
    }
}

void partial_order::copy(const partial_order& order, 
		         map< indicator_term *, indicator_term * > old2new)
{
    map<indicator_term *, vector<indicator_term * > >::const_iterator i;
    map<indicator_term *, int >::const_iterator j;

    for (j = order.pred.begin(); j != order.pred.end(); ++j)
	pred[old2new[(*j).first]] = (*j).second;

    for (i = order.lt.begin(); i != order.lt.end(); ++i) {
	for (int j = 0; j < (*i).second.size(); ++j)
	    lt[old2new[(*i).first]].push_back(old2new[(*i).second[j]]);
    }
    for (i = order.le.begin(); i != order.le.end(); ++i) {
	for (int j = 0; j < (*i).second.size(); ++j)
	    le[old2new[(*i).first]].push_back(old2new[(*i).second[j]]);
    }
    for (i = order.eq.begin(); i != order.eq.end(); ++i) {
	for (int j = 0; j < (*i).second.size(); ++j)
	    eq[old2new[(*i).first]].push_back(old2new[(*i).second[j]]);
    }
    for (i = order.pending.begin(); i != order.pending.end(); ++i) {
	for (int j = 0; j < (*i).second.size(); ++j)
	    pending[old2new[(*i).first]].push_back(old2new[(*i).second[j]]);
    }
}

static void interval_minmax(Polyhedron *I, int *min, int *max)
{
    assert(I->Dimension == 1);
    *min = 1;
    *max = -1;
    POL_ENSURE_VERTICES(I);
    for (int i = 0; i < I->NbRays; ++i) {
	if (value_pos_p(I->Ray[i][1]))
	    *max = 1;
	else if (value_neg_p(I->Ray[i][1]))
	    *min = -1;
	else {
	    if (*max < 0)
		*max = 0;
	    if (*min > 0)
		*min = 0;
	}
    }
}

static void interval_minmax(Polyhedron *D, Matrix *T, int *min, int *max, 
			    unsigned MaxRays)
{
    Polyhedron *I = Polyhedron_Image(D, T, MaxRays);
    if (MaxRays & POL_INTEGER)
	I = DomainConstraintSimplify(I, MaxRays);
    if (emptyQ2(I)) {
	Polyhedron_Free(I);
	I = Polyhedron_Image(D, T, MaxRays);
    }
    interval_minmax(I, min, max);
    Polyhedron_Free(I);
}

struct max_term {
    unsigned dim;
    EDomain *domain;
    vector<evalue *> max;

    void print(ostream& os, char **p) const;
    void resolve_existential_vars() const;
    void substitute(Matrix *T, unsigned MaxRays);
    Vector *eval(Value *val, unsigned MaxRays) const;

    ~max_term() {
	for (int i = 0; i < max.size(); ++i) {
	    free_evalue_refs(max[i]);
	    delete max[i];
	}
	delete domain;
    }
};

/*
 * Project on first dim dimensions
 */
Polyhedron* Polyhedron_Project_Initial(Polyhedron *P, int dim)
{
    int i;
    Matrix *T;
    Polyhedron *I;

    if (P->Dimension == dim)
	return Polyhedron_Copy(P);

    T = Matrix_Alloc(dim+1, P->Dimension+1);
    for (i = 0; i < dim; ++i)
	value_set_si(T->p[i][i], 1);
    value_set_si(T->p[dim][P->Dimension], 1);
    I = Polyhedron_Image(P, T, P->NbConstraints);
    Matrix_Free(T);
    return I;
}

struct indicator {
    vector<indicator_term*> term;
    indicator_constructor& ic;
    partial_order order;
    EDomain	 *D;
    Polyhedron	 *P;
    Param_Domain *PD;
    barvinok_options	*options;

    indicator(indicator_constructor& ic, Param_Domain *PD, EDomain *D,
	      barvinok_options *options) :
	ic(ic), PD(PD), D(D), order(this), options(options), P(NULL) {}
    indicator(const indicator& ind, EDomain *D) :
	    ic(ind.ic), PD(ind.PD), D(NULL), order(this), options(ind.options),
	    P(Polyhedron_Copy(ind.P)) {
	map< indicator_term *, indicator_term * > old2new;
	for (int i = 0; i < ind.term.size(); ++i) {
	    indicator_term *it = new indicator_term(*ind.term[i]);
	    old2new[ind.term[i]] = it;
	    term.push_back(it);
	}
	order.copy(ind.order, old2new);
	set_domain(D);
    }
    ~indicator() {
	for (int i = 0; i < term.size(); ++i)
	    delete term[i];
	if (D)
	    delete D;
	if (P)
	    Polyhedron_Free(P);
    }

    void set_domain(EDomain *D) {
	if (this->D)
	    delete this->D;
	this->D = D;
	int nparam = ic.P->Dimension - ic.vertex.length();
	Polyhedron *Q = Polyhedron_Project_Initial(D->D, nparam);
	Q = DomainConstraintSimplify(Q, options->MaxRays);
	if (!P || !PolyhedronIncludes(Q, P))
	    reduce_in_domain(Q);
	if (P)
	    Polyhedron_Free(P);
	P = Q;
    }

    void add(const indicator_term* it);
    void remove(indicator_term* it);
    void remove_initial_rational_vertices();
    void expand_rational_vertex(indicator_term *initial);

    void print(ostream& os, char **p);
    void simplify();
    void peel(int i, int j);
    void combine(indicator_term *a, indicator_term *b);
    void substitute(evalue *equation);
    void reduce_in_domain(Polyhedron *D);
    bool handle_equal_numerators(indicator_term *base);

    max_term* create_max_term(indicator_term *it);
};

max_term* indicator::create_max_term(indicator_term *it)
{
    int dim = it->den.NumCols();
    int nparam = ic.P->Dimension - ic.vertex.length();
    max_term *maximum = new max_term;
    maximum->dim = nparam;
    maximum->domain = new EDomain(D);
    for (int j = 0; j < dim; ++j) {
	evalue *E = new evalue;
	value_init(E->d);
	evalue_copy(E, it->vertex[j]);
	if (evalue_frac2floor_in_domain(E, D->D))
	    reduce_evalue(E);
	maximum->max.push_back(E);
    }
    return maximum;
}

/* compute a - b */
static evalue *ediff(const evalue *a, const evalue *b)
{
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
    evalue *diff = new evalue;
    value_init(diff->d);
    evalue_copy(diff, b);
    emul(&mone, diff);
    eadd(a, diff);
    reduce_evalue(diff);
    free_evalue_refs(&mone); 
    return diff;
}

static order_sign evalue_sign(evalue *diff, EDomain *D, unsigned MaxRays)
{
    order_sign sign = order_eq;
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
    int len = 1 + D->D->Dimension + 1;
    Vector *c = Vector_Alloc(len);
    Matrix *T = Matrix_Alloc(2, len-1);

    int fract = evalue2constraint(D, diff, c->p, len);
    Vector_Copy(c->p+1, T->p[0], len-1);
    value_assign(T->p[1][len-2], c->p[0]);

    int min, max;
    interval_minmax(D->D, T, &min, &max, MaxRays);
    if (max < 0)
	sign = order_lt;
    else {
	if (fract) {
	    emul(&mone, diff);
	    evalue2constraint(D, diff, c->p, len);
	    emul(&mone, diff);
	    Vector_Copy(c->p+1, T->p[0], len-1);
	    value_assign(T->p[1][len-2], c->p[0]);

	    int negmin, negmax;
	    interval_minmax(D->D, T, &negmin, &negmax, MaxRays);
	    min = -negmax;
	}
	if (min > 0)
	    sign = order_gt;
	else if (max == 0 && min == 0)
	    sign = order_eq;
	else if (min < 0 && max > 0)
	    sign = order_unknown;
	else if (min < 0)
	    sign = order_le;
	else
	    sign = order_ge;
    }

    Matrix_Free(T);
    Vector_Free(c);
    free_evalue_refs(&mone); 

    return sign;
}

order_sign partial_order::compare(indicator_term *a, indicator_term *b)
{
    unsigned dim = a->den.NumCols();
    order_sign sign = order_eq;
    EDomain *D = ind->D;
    unsigned MaxRays = ind->options->MaxRays;
    if (MaxRays & POL_INTEGER && (a->sign == 0 || b->sign == 0))
	MaxRays = 0;

    for (int k = 0; k < dim; ++k) {
	/* compute a->vertex[k] - b->vertex[k] */
	evalue *diff = ediff(a->vertex[k], b->vertex[k]);
	order_sign diff_sign = evalue_sign(diff, D, MaxRays);

	if (diff_sign == order_lt) {
	    if (sign == order_eq || sign == order_le)
		sign = order_lt;
	    else
		sign = order_unknown;
	    free_evalue_refs(diff); 
	    delete diff;
	    break;
	}
	if (diff_sign == order_gt) {
	    if (sign == order_eq || sign == order_ge)
		sign = order_gt;
	    else
		sign = order_unknown;
	    free_evalue_refs(diff); 
	    delete diff;
	    break;
	}
	if (diff_sign == order_eq) {
	    if (D == ind->D && !EVALUE_IS_ZERO(*diff))
		ind->substitute(diff);
	    free_evalue_refs(diff); 
	    delete diff;
	    continue;
	}
	if ((diff_sign == order_unknown) ||
	    ((diff_sign == order_lt || diff_sign == order_le) && sign == order_ge) ||
	    ((diff_sign == order_gt || diff_sign == order_ge) && sign == order_le)) {
	    free_evalue_refs(diff); 
	    delete diff;
	    sign = order_unknown;
	    break;
	}

	sign = diff_sign;

	Matrix *M;
	vector<EDomain_floor *> new_floors;
	M = D->add_ge_constraint(diff, new_floors);
	value_set_si(M->p[M->NbRows-1][0], 0);
	Polyhedron *D2 = Constraints2Polyhedron(M, MaxRays);
	EDomain *EDeq = new EDomain(D2, D, new_floors);
	Polyhedron_Free(D2);
	Matrix_Free(M);
	for (int i = 0; i < new_floors.size(); ++i)
	    EDomain_floor::unref(new_floors[i]);

	if (D != ind->D)
	    delete D;
	D = EDeq;

	free_evalue_refs(diff); 
	delete diff;
    }

    if (D != ind->D)
	delete D;

    return sign;
}

bool partial_order::compared(indicator_term* a, indicator_term* b)
{
    map<indicator_term *, vector<indicator_term * > >::iterator j;

    j = lt.find(a);
    if (j != lt.end() && find(lt[a].begin(), lt[a].end(), b) != lt[a].end())
	return true;

    j = le.find(a);
    if (j != le.end() && find(le[a].begin(), le[a].end(), b) != le[a].end())
	return true;

    return false;
}

void partial_order::add(indicator_term* it, std::set<indicator_term *> *filter)
{
    if (eq.find(it) != eq.end() && eq[it].size() == 1)
	return;

    int it_pred = filter ? pred[it] : 0;

    map<indicator_term *, int >::iterator i;
    for (i = pred.begin(); i != pred.end(); ++i) {
	if ((*i).second != 0)
	    continue;
	if (eq.find((*i).first) != eq.end() && eq[(*i).first].size() == 1)
	    continue;
	if (filter) {
	    if ((*i).first == it)
		continue;
	    if (filter->find((*i).first) == filter->end())
		continue;
	    if (compared((*i).first, it))
		continue;
	}
	order_sign sign = compare(it, (*i).first);
	if (sign == order_lt) {
	    lt[it].push_back((*i).first);
	    (*i).second++;
	} else if (sign == order_le) {
	    le[it].push_back((*i).first);
	    (*i).second++;
	} else if (sign == order_eq) {
	    pred[it] = it_pred;
	    set_equal(it, (*i).first);
	    return;
	} else if (sign == order_gt) {
	    pending[(*i).first].push_back(it);
	    lt[(*i).first].push_back(it);
	    ++it_pred;
	} else if (sign == order_ge) {
	    pending[(*i).first].push_back(it);
	    le[(*i).first].push_back(it);
	    ++it_pred;
	}
    }
    pred[it] = it_pred;
}

void partial_order::remove(indicator_term* it)
{
    std::set<indicator_term *> filter;
    map<indicator_term *, vector<indicator_term * > >::iterator i;

    assert(pred[it] == 0);

    i = eq.find(it);
    if (i != eq.end()) {
	assert(eq[it].size() >= 1);
	indicator_term *base;
	if (eq[it].size() == 1) {
	    base = eq[it][0];
	    eq.erase(i);

	    vector<indicator_term * >::iterator j;
	    j = find(eq[base].begin(), eq[base].end(), it);
	    assert(j != eq[base].end());
	    eq[base].erase(j);
	} else {
	    /* "it" may no longer be the smallest, since the order
	     * structure may have been copied from another one.
	     */
	    sort(eq[it].begin()+1, eq[it].end());
	    assert(eq[it][0] == it);
	    eq[it].erase(eq[it].begin());
	    base = eq[it][0];
	    eq[base] = eq[it];
	    eq.erase(i);

	    for (int j = 1; j < eq[base].size(); ++j)
		eq[eq[base][j]][0] = base;

	    i = lt.find(it);
	    if (i != lt.end()) {
		lt[base] = lt[it];
		lt.erase(i);
	    }
	
	    i = le.find(it);
	    if (i != le.end()) {
		le[base] = le[it];
		le.erase(i);
	    }
	
	    i = pending.find(it);
	    if (i != pending.end()) {
		pending[base] = pending[it];
		pending.erase(i);
	    }
	}

	if (eq[base].size() == 1)
	    eq.erase(base);

	map<indicator_term *, int >::iterator j;
	j = pred.find(it);
	pred.erase(j);

	return;
    }

    i = lt.find(it);
    if (i != lt.end()) {
	for (int j = 0; j < lt[it].size(); ++j) {
	    filter.insert(lt[it][j]);
	    pred[lt[it][j]]--;
	}
	lt.erase(i);
    }

    i = le.find(it);
    if (i != le.end()) {
	for (int j = 0; j < le[it].size(); ++j) {
	    filter.insert(le[it][j]);
	    pred[le[it][j]]--;
	}
	le.erase(i);
    }

    map<indicator_term *, int >::iterator j;
    j = pred.find(it);
    pred.erase(j);

    i = pending.find(it);
    if (i != pending.end()) {
	for (int j = 0; j < pending[it].size(); ++j) {
	    filter.erase(pending[it][j]);
	    add(pending[it][j], &filter);
	}
	pending.erase(i);
    }
}

void partial_order::print(ostream& os, char **p)
{
    map<indicator_term *, vector<indicator_term * > >::iterator i;
    for (i = lt.begin(); i != lt.end(); ++i) {
	(*i).first->print(os, p);
	assert(pred.find((*i).first) != pred.end());
	os << "(" << pred[(*i).first] << ")";
	os << " < ";
	for (int j = 0; j < (*i).second.size(); ++j) {
	    if (j)
		os << ", ";
	    (*i).second[j]->print(os, p);
	    assert(pred.find((*i).second[j]) != pred.end());
	    os << "(" << pred[(*i).second[j]] << ")";
	}
	os << endl;
    }
    for (i = le.begin(); i != le.end(); ++i) {
	(*i).first->print(os, p);
	assert(pred.find((*i).first) != pred.end());
	os << "(" << pred[(*i).first] << ")";
	os << " <= ";
	for (int j = 0; j < (*i).second.size(); ++j) {
	    if (j)
		os << ", ";
	    (*i).second[j]->print(os, p);
	    assert(pred.find((*i).second[j]) != pred.end());
	    os << "(" << pred[(*i).second[j]] << ")";
	}
	os << endl;
    }
    for (i = eq.begin(); i != eq.end(); ++i) {
	if ((*i).second.size() <= 1)
	    continue;
	(*i).first->print(os, p);
	assert(pred.find((*i).first) != pred.end());
	os << "(" << pred[(*i).first] << ")";
	for (int j = 1; j < (*i).second.size(); ++j) {
	    if (j)
		os << " = ";
	    (*i).second[j]->print(os, p);
	    assert(pred.find((*i).second[j]) != pred.end());
	    os << "(" << pred[(*i).second[j]] << ")";
	}
	os << endl;
    }
    for (i = pending.begin(); i != pending.end(); ++i) {
	os << "pending on ";
	(*i).first->print(os, p);
	assert(pred.find((*i).first) != pred.end());
	os << "(" << pred[(*i).first] << ")";
	os << ": ";
	for (int j = 0; j < (*i).second.size(); ++j) {
	    if (j)
		os << ", ";
	    (*i).second[j]->print(os, p);
	    assert(pred.find((*i).second[j]) != pred.end());
	    os << "(" << pred[(*i).second[j]] << ")";
	}
	os << endl;
    }
}

void indicator::add(const indicator_term* it)
{
    indicator_term *nt = new indicator_term(*it);
    nt->reduce_in_domain(P ? P : D->D);
    term.push_back(nt);
    order.add(nt, NULL);
    assert(term.size() == order.pred.size());
}

void indicator::remove(indicator_term* it)
{
    vector<indicator_term *>::iterator i;
    i = find(term.begin(), term.end(), it);
    assert(i!= term.end());
    order.remove(it);
    term.erase(i);
    assert(term.size() == order.pred.size());
    delete it;
}

void indicator::expand_rational_vertex(indicator_term *initial)
{
    int pos = initial->pos;
    remove(initial);
    if (ic.terms[pos].size() == 0) {
	Param_Vertices *V;
	FORALL_PVertex_in_ParamPolyhedron(V, PD, ic.PP) // _i is internal counter
	    if (_i == pos) {
		ic.decompose_at_vertex(V, pos, options->MaxRays);
		break;
	    }
	END_FORALL_PVertex_in_ParamPolyhedron;
    }
    for (int j = 0; j < ic.terms[pos].size(); ++j)
	add(ic.terms[pos][j]);
}

void indicator::remove_initial_rational_vertices()
{
    do {
	indicator_term *initial = NULL;
	map<indicator_term *, int >::iterator i;
	for (i = order.pred.begin(); i != order.pred.end(); ++i) {
	    if ((*i).second != 0)
		continue;
	    if ((*i).first->sign != 0)
		continue;
	    if (order.eq.find((*i).first) != order.eq.end() &&
		order.eq[(*i).first].size() <= 1)
		continue;
	    initial = (*i).first;
	    break;
	}
	if (!initial)
	    break;
	expand_rational_vertex(initial);
    } while(1);
}

void indicator::reduce_in_domain(Polyhedron *D)
{
    for (int i = 0; i < term.size(); ++i)
	term[i]->reduce_in_domain(D);
}

void indicator::print(ostream& os, char **p)
{
    assert(term.size() == order.pred.size());
    for (int i = 0; i < term.size(); ++i) {
	term[i]->print(os, p);
	os << endl;
    }
    order.print(os, p);
}

/* Remove pairs of opposite terms */
void indicator::simplify()
{
    for (int i = 0; i < term.size(); ++i) {
	for (int j = i+1; j < term.size(); ++j) {
	    if (term[i]->sign + term[j]->sign != 0)
		continue;
	    if (term[i]->den != term[j]->den)
		continue;
	    int k;
	    for (k = 0; k < term[i]->den.NumCols(); ++k)
		if (!eequal(term[i]->vertex[k], term[j]->vertex[k]))
		    break;
	    if (k < term[i]->den.NumCols())
		continue;
	    delete term[i];
	    delete term[j];
	    term.erase(term.begin()+j);
	    term.erase(term.begin()+i);
	    --i;
	    break;
	}
    }
}

void indicator::peel(int i, int j)
{
    if (j < i) {
	int tmp = i;
	i = j;
	j = tmp;
    }
    assert (i < j);
    int dim = term[i]->den.NumCols();

    mat_ZZ common;
    mat_ZZ rest_i;
    mat_ZZ rest_j;
    int n_common = 0, n_i = 0, n_j = 0;

    common.SetDims(min(term[i]->den.NumRows(), term[j]->den.NumRows()), dim);
    rest_i.SetDims(term[i]->den.NumRows(), dim);
    rest_j.SetDims(term[j]->den.NumRows(), dim);

    int k, l;
    for (k = 0, l = 0; k < term[i]->den.NumRows() && l < term[j]->den.NumRows(); ) {
	int s = lex_cmp(term[i]->den[k], term[j]->den[l]);
	if (s == 0) {
	    common[n_common++] = term[i]->den[k];
	    ++k;
	    ++l;
	} else if (s < 0)
	    rest_i[n_i++] = term[i]->den[k++];
	else
	    rest_j[n_j++] = term[j]->den[l++];
    }
    while (k < term[i]->den.NumRows())
	rest_i[n_i++] = term[i]->den[k++];
    while (l < term[j]->den.NumRows())
	rest_j[n_j++] = term[j]->den[l++];
    common.SetDims(n_common, dim);
    rest_i.SetDims(n_i, dim);
    rest_j.SetDims(n_j, dim);

    for (k = 0; k <= n_i; ++k) {
	indicator_term *it = new indicator_term(*term[i]);
	it->den.SetDims(n_common + k, dim);
	for (l = 0; l < n_common; ++l)
	    it->den[l] = common[l];
	for (l = 0; l < k; ++l)
	    it->den[n_common+l] = rest_i[l];
	lex_order_rows(it->den);
	if (k)
	    for (l = 0; l < dim; ++l)
		evalue_add_constant(it->vertex[l], rest_i[k-1][l]);
	term.push_back(it);
    }

    for (k = 0; k <= n_j; ++k) {
	indicator_term *it = new indicator_term(*term[j]);
	it->den.SetDims(n_common + k, dim);
	for (l = 0; l < n_common; ++l)
	    it->den[l] = common[l];
	for (l = 0; l < k; ++l)
	    it->den[n_common+l] = rest_j[l];
	lex_order_rows(it->den);
	if (k)
	    for (l = 0; l < dim; ++l)
		evalue_add_constant(it->vertex[l], rest_j[k-1][l]);
	term.push_back(it);
    }
    term.erase(term.begin()+j);
    term.erase(term.begin()+i);
}

void indicator::combine(indicator_term *a, indicator_term *b)
{
    int dim = a->den.NumCols();

    mat_ZZ common;
    mat_ZZ rest_i;
    mat_ZZ rest_j;
    int n_common = 0, n_i = 0, n_j = 0;

    common.SetDims(min(a->den.NumRows(), b->den.NumRows()), dim);
    rest_i.SetDims(a->den.NumRows(), dim);
    rest_j.SetDims(b->den.NumRows(), dim);

    int k, l;
    for (k = 0, l = 0; k < a->den.NumRows() && l < b->den.NumRows(); ) {
	int s = lex_cmp(a->den[k], b->den[l]);
	if (s == 0) {
	    common[n_common++] = a->den[k];
	    ++k;
	    ++l;
	} else if (s < 0)
	    rest_i[n_i++] = a->den[k++];
	else
	    rest_j[n_j++] = b->den[l++];
    }
    while (k < a->den.NumRows())
	rest_i[n_i++] = a->den[k++];
    while (l < b->den.NumRows())
	rest_j[n_j++] = b->den[l++];
    common.SetDims(n_common, dim);
    rest_i.SetDims(n_i, dim);
    rest_j.SetDims(n_j, dim);

    assert(n_i < 30);
    assert(n_j < 30);

    for (k = (1 << n_i)-1; k >= 0; --k) {
	indicator_term *it = k ? new indicator_term(*b) : b;
	it->den.SetDims(n_common + n_i + n_j, dim);
	for (l = 0; l < n_common; ++l)
	    it->den[l] = common[l];
	for (l = 0; l < n_i; ++l)
	    it->den[n_common+l] = rest_i[l];
	for (l = 0; l < n_j; ++l)
	    it->den[n_common+n_i+l] = rest_j[l];
	lex_order_rows(it->den);
	int change = 0;
	for (l = 0; l < n_i; ++l) {
	    if (!(k & (1 <<l)))
		continue;
	    change ^= 1;
	    for (int m = 0; m < dim; ++m)
		evalue_add_constant(it->vertex[m], rest_i[l][m]);
	}
	if (change)
	    it->sign = -it->sign;
	if (it != b) {
	    term.push_back(it);
	    order.add(it, NULL);
	}
    }

    for (k = (1 << n_j)-1; k >= 0; --k) {
	indicator_term *it = k ? new indicator_term(*a) : a;
	it->den.SetDims(n_common + n_i + n_j, dim);
	for (l = 0; l < n_common; ++l)
	    it->den[l] = common[l];
	for (l = 0; l < n_i; ++l)
	    it->den[n_common+l] = rest_i[l];
	for (l = 0; l < n_j; ++l)
	    it->den[n_common+n_i+l] = rest_j[l];
	lex_order_rows(it->den);
	int change = 0;
	for (l = 0; l < n_j; ++l) {
	    if (!(k & (1 <<l)))
		continue;
	    change ^= 1;
	    for (int m = 0; m < dim; ++m)
		evalue_add_constant(it->vertex[m], rest_j[l][m]);
	}
	if (change)
	    it->sign = -it->sign;
	if (it != a) {
	    term.push_back(it);
	    order.add(it, NULL);
	}
    }
}

bool indicator::handle_equal_numerators(indicator_term *base)
{
    for (int i = 0; i < order.eq[base].size(); ++i) {
	for (int j = i+1; j < order.eq[base].size(); ++j) {
	    if (order.eq[base][i]->is_opposite(order.eq[base][j])) {
		remove(order.eq[base][j]);
		remove(i ? order.eq[base][i] : base);
		return true;
	    }
	}
    }
    for (int j = 1; j < order.eq[base].size(); ++j)
	if (order.eq[base][j]->sign != base->sign) {
	    combine(base, order.eq[base][j]);
	    return true;
	}
    return false;
}

void indicator::substitute(evalue *equation)
{
    evalue *fract = NULL;
    evalue *val = new evalue;
    value_init(val->d);
    evalue_copy(val, equation);

    evalue factor;
    value_init(factor.d);
    value_init(factor.x.n);

    evalue *e;
    for (e = val; value_zero_p(e->d) && e->x.p->type != fractional;
	 e = &e->x.p->arr[type_offset(e->x.p)])
	;

    if (value_zero_p(e->d) && e->x.p->type == fractional)
	fract = &e->x.p->arr[0];
    else {
	for (e = val; value_zero_p(e->d) && e->x.p->type != polynomial;
	     e = &e->x.p->arr[type_offset(e->x.p)])
	    ;
	assert(value_zero_p(e->d) && e->x.p->type == polynomial);
    }

    int offset = type_offset(e->x.p);

    assert(value_notzero_p(e->x.p->arr[offset+1].d));
    assert(value_notzero_p(e->x.p->arr[offset+1].x.n));
    if (value_neg_p(e->x.p->arr[offset+1].x.n)) {
	value_oppose(factor.d, e->x.p->arr[offset+1].x.n);
	value_assign(factor.x.n, e->x.p->arr[offset+1].d);
    } else {
	value_assign(factor.d, e->x.p->arr[offset+1].x.n);
	value_oppose(factor.x.n, e->x.p->arr[offset+1].d);
    }

    free_evalue_refs(&e->x.p->arr[offset+1]);
    enode *p = e->x.p;
    value_clear(e->d);
    *e = e->x.p->arr[offset];

    emul(&factor, val);

    if (fract) {
	for (int i = 0; i < term.size(); ++i)
	    term[i]->substitute(fract, val);

	free_evalue_refs(&p->arr[0]);
    } else {
	for (int i = 0; i < term.size(); ++i)
	    term[i]->substitute(p->pos, val);
    }

    free_evalue_refs(&factor);
    free_evalue_refs(val);
    delete val;

    free(p);
}

static void print_varlist(ostream& os, int n, char **names)
{
    int i;
    os << "[";
    for (i = 0; i < n; ++i) {
	if (i)
	    os << ",";
	os << names[i];
    }
    os << "]";
}

static void print_term(ostream& os, Value v, int pos, int dim,
		        char **names, int *first)
{
    if (value_zero_p(v)) {
	if (first && *first && pos >= dim)
	    os << "0";
	return;
    }

    if (first) {
	if (!*first && value_pos_p(v))
	    os << "+";
	*first = 0;
    }
    if (pos < dim) {
	if (value_mone_p(v)) {
	    os << "-";
	} else if (!value_one_p(v))
	    os << VALUE_TO_INT(v);
	os << names[pos];
    } else
	os << VALUE_TO_INT(v);
}

/* We put all possible existentially quantified variables at the back
 * and so if any equalities exist between these variables and the
 * other variables, then PolyLib will replace all occurrences of some
 * of the other variables by some existentially quantified variables.
 * We want the output to have as few as possible references to the
 * existentially quantified variables, so we undo what PolyLib did here.
 */
void resolve_existential_vars(Polyhedron *domain, unsigned dim)
{
    int last = domain->NbEq - 1;
    /* Matrix "view" of domain for ExchangeRows */
    Matrix M;
    M.NbRows = domain->NbConstraints;
    M.NbColumns = domain->Dimension+2;
    M.p_Init = domain->p_Init;
    M.p = domain->Constraint;
    Value mone;
    value_init(mone);
    value_set_si(mone, -1);
    for (int e = domain->Dimension-1; e >= dim; --e) {
	int r;
	for (r = last; r >= 0; --r)
	    if (value_notzero_p(domain->Constraint[r][1+e]))
		break;
	if (r < 0)
	    continue;
	if (r != last)
	    ExchangeRows(&M, r, last);

	/* Combine uses the coefficient as a multiplier, so it must
	 * be positive, since we are modifying an inequality.
	 */
	if (value_neg_p(domain->Constraint[last][1+e]))
	    Vector_Scale(domain->Constraint[last]+1, domain->Constraint[last]+1,
			 mone, domain->Dimension+1);

	for (int c = 0; c < domain->NbConstraints; ++c) {
	    if (c == last)
		continue;
	    if (value_notzero_p(domain->Constraint[c][1+e]))
		Combine(domain->Constraint[c], domain->Constraint[last],
			domain->Constraint[c], 1+e, domain->Dimension+1);
	}
	--last;
    }
    value_clear(mone);
}

void max_term::resolve_existential_vars() const
{
    ::resolve_existential_vars(domain->D, dim);
}

void max_term::print(ostream& os, char **p) const
{
    char **names = p;
    if (dim < domain->D->Dimension) {
	resolve_existential_vars();
	names = new char * [domain->D->Dimension];
	int i;
	for (i = 0; i < dim; ++i)
	    names[i] = p[i];
	for ( ; i < domain->D->Dimension; ++i) {
	    names[i] = new char[10];
	    snprintf(names[i], 10, "a%d", i - dim);
	}
    }

    Value tmp;
    value_init(tmp);
    os << "{ ";
    print_varlist(os, dim, p);
    os << " -> ";
    os << "[";
    for (int i = 0; i < max.size(); ++i) {
	if (i)
	    os << ",";
	evalue_print(os, max[i], p);
    }
    os << "]";
    os << " : ";
    if (dim < domain->D->Dimension) {
	os << "exists ";
	print_varlist(os, domain->D->Dimension-dim, names+dim);
	os << " : ";
    }
    for (int i = 0; i < domain->D->NbConstraints; ++i) {
	int first = 1;
	int v = First_Non_Zero(domain->D->Constraint[i]+1, domain->D->Dimension);
	if (v == -1)
	    continue;
	if (i)
	    os << " && ";
	if (value_pos_p(domain->D->Constraint[i][v+1])) {
	    print_term(os, domain->D->Constraint[i][v+1], v, domain->D->Dimension,
		       names, NULL);
	    if (value_zero_p(domain->D->Constraint[i][0]))
		os << " = ";
	    else
		os << " >= ";
	    for (int j = v+1; j <= domain->D->Dimension; ++j) {
		value_oppose(tmp, domain->D->Constraint[i][1+j]);
		print_term(os, tmp, j, domain->D->Dimension,
			   names, &first);
	    }
	} else {
	    value_oppose(tmp, domain->D->Constraint[i][1+v]);
	    print_term(os, tmp, v, domain->D->Dimension,
		       names, NULL);
	    if (value_zero_p(domain->D->Constraint[i][0]))
		os << " = ";
	    else
		os << " <= ";
	    for (int j = v+1; j <= domain->D->Dimension; ++j)
		print_term(os, domain->D->Constraint[i][1+j], j, domain->D->Dimension,
			   names, &first);
	}
    }
    os << " }" << endl;
    value_clear(tmp);

    if (dim < domain->D->Dimension) {
	for (int i = dim; i < domain->D->Dimension; ++i)
	    delete [] names[i];
	delete [] names;
    }
}

static void evalue_substitute(evalue *e, evalue **subs)
{
    evalue *v;

    if (value_notzero_p(e->d))
	return;

    enode *p = e->x.p;
    for (int i = 0; i < p->size; ++i)
	evalue_substitute(&p->arr[i], subs);

    if (p->type == polynomial)
	v = subs[p->pos-1];
    else {
	v = new evalue;
	value_init(v->d);
	value_set_si(v->d, 0);
	v->x.p = new_enode(p->type, 3, -1);
	value_clear(v->x.p->arr[0].d);
	v->x.p->arr[0] = p->arr[0];
	evalue_set_si(&v->x.p->arr[1], 0, 1);
	evalue_set_si(&v->x.p->arr[2], 1, 1);
    }

    int offset = type_offset(p);

    for (int i = p->size-1; i >= offset+1; i--) {
	emul(v, &p->arr[i]);
	eadd(&p->arr[i], &p->arr[i-1]);
	free_evalue_refs(&(p->arr[i]));
    }

    if (p->type != polynomial) {
	free_evalue_refs(v);
	delete v;
    }

    value_clear(e->d);
    *e = p->arr[offset];
    free(p);
}

/* "align" matrix to have nrows by inserting
 * the necessary number of rows and an equal number of columns at the end
 * right before the constant row/column
 */
static Matrix *align_matrix_initial(Matrix *M, int nrows)
{
    int i;
    int newrows = nrows - M->NbRows;
    Matrix *M2 = Matrix_Alloc(nrows, newrows + M->NbColumns);
    for (i = 0; i < newrows; ++i)
	value_set_si(M2->p[M->NbRows-1+i][M->NbColumns-1+i], 1);
    for (i = 0; i < M->NbRows-1; ++i) {
	Vector_Copy(M->p[i], M2->p[i], M->NbColumns-1);
	value_assign(M2->p[i][M2->NbColumns-1], M->p[i][M->NbColumns-1]);
    }
    value_assign(M2->p[M2->NbRows-1][M2->NbColumns-1],
		 M->p[M->NbRows-1][M->NbColumns-1]);
    return M2;
}

/* T maps the compressed parameters to the original parameters,
 * while this max_term is based on the compressed parameters
 * and we want get the original parameters back.
 */
void max_term::substitute(Matrix *T, unsigned MaxRays)
{
    int nexist = domain->D->Dimension - (T->NbColumns-1);
    Matrix *M = align_matrix_initial(T, T->NbRows+nexist);

    Polyhedron *D = DomainImage(domain->D, M, MaxRays);
    Polyhedron_Free(domain->D);
    domain->D = D;
    Matrix_Free(M);

    assert(T->NbRows == T->NbColumns);
    Matrix *T2 = Matrix_Copy(T);
    Matrix *inv = Matrix_Alloc(T->NbColumns, T->NbRows);
    int ok = Matrix_Inverse(T2, inv);
    Matrix_Free(T2);
    assert(ok);

    evalue denom;
    value_init(denom.d);
    value_init(denom.x.n);
    value_set_si(denom.x.n, 1);
    value_assign(denom.d, inv->p[inv->NbRows-1][inv->NbColumns-1]);

    vec_ZZ v;
    v.SetLength(inv->NbColumns);
    evalue* subs[inv->NbRows-1];
    for (int i = 0; i < inv->NbRows-1; ++i) {
	values2zz(inv->p[i], v, v.length());
	subs[i] = multi_monom(v);
	emul(&denom, subs[i]);
    }
    free_evalue_refs(&denom);

    for (int i = 0; i < max.size(); ++i) {
	evalue_substitute(max[i], subs);
	reduce_evalue(max[i]);
    }

    for (int i = 0; i < inv->NbRows-1; ++i) {
	free_evalue_refs(subs[i]);
	delete subs[i];
    }
    Matrix_Free(inv);
}

int Last_Non_Zero(Value *p, unsigned len)
{
    for (int i = len-1; i >= 0; --i)
	if (value_notzero_p(p[i]))
	    return i;
    return -1;
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

bool in_domain(Polyhedron *P, Value *val, unsigned dim, unsigned MaxRays)
{
    int nexist = P->Dimension - dim;
    int last[P->NbConstraints];
    Value tmp, min, max;
    Vector *all_val = Vector_Alloc(P->Dimension+1);
    bool in = false;
    int i;
    int alternate;

    resolve_existential_vars(P, dim);

    Vector_Copy(val, all_val->p, dim);
    value_set_si(all_val->p[P->Dimension], 1);

    value_init(tmp);
    for (int i = 0; i < P->NbConstraints; ++i) {
	last[i] = Last_Non_Zero(P->Constraint[i]+1+dim, nexist);
	if (last[i] == -1) {
	    Inner_Product(P->Constraint[i]+1, all_val->p, P->Dimension+1, &tmp);
	    if (value_neg_p(tmp))
		goto out;
	    if (i < P->NbEq && value_pos_p(tmp))
		goto out;
	}
    }

    value_init(min);
    value_init(max);
    alternate = nexist - 1;
    for (i = 0; i < nexist; ++i) {
	bool min_set = false;
	bool max_set = false;
	for (int j = 0; j < P->NbConstraints; ++j) {
	    if (last[j] != i)
		continue;
	    Inner_Product(P->Constraint[j]+1, all_val->p, P->Dimension+1, &tmp);
	    value_oppose(tmp, tmp);
	    if (j < P->NbEq) {
		if (!mpz_divisible_p(tmp, P->Constraint[j][1+dim+i]))
		    goto out2;
		value_division(tmp, tmp, P->Constraint[j][1+dim+i]);
		if (!max_set || value_lt(tmp, max)) {
		    max_set = true;
		    value_assign(max, tmp);
		}
		if (!min_set || value_gt(tmp, min)) {
		    min_set = true;
		    value_assign(min, tmp);
		}
	    } else {
		if (value_pos_p(P->Constraint[j][1+dim+i])) {
		    mpz_cdiv_q(tmp, tmp, P->Constraint[j][1+dim+i]);
		    if (!min_set || value_gt(tmp, min)) {
			min_set = true;
			value_assign(min, tmp);
		    }
		} else {
		    mpz_fdiv_q(tmp, tmp, P->Constraint[j][1+dim+i]);
		    if (!max_set || value_lt(tmp, max)) {
			max_set = true;
			value_assign(max, tmp);
		    }
		}
	    }
	}
	/* Move another existential variable in current position */
	if (!max_set || !min_set) {
	    if (!(alternate > i)) {
		Matrix *M = Matrix_Alloc(dim+i, 1+P->Dimension+1);
		for (int j = 0; j < dim+i; ++j) {
		    value_set_si(M->p[j][1+j], -1);
		    value_assign(M->p[j][1+P->Dimension], all_val->p[j]);
		}
		Polyhedron *Q = AddConstraints(M->p[0], dim+i, P, MaxRays);
		Matrix_Free(M);
		Q = DomainConstraintSimplify(Q, MaxRays);
		Vector *sample = Polyhedron_Sample(Q, MaxRays);
		in = !!sample;
		if (sample)
		    Vector_Free(sample);
		Polyhedron_Free(Q);
		goto out2;
	    }
	    assert(alternate > i);
	    SwapColumns(P, 1+dim+i, 1+dim+alternate);
	    resolve_existential_vars(P, dim);
	    for (int j = 0; j < P->NbConstraints; ++j) {
		if (j >= P->NbEq && last[j] < i)
		    continue;
		last[j] = Last_Non_Zero(P->Constraint[j]+1+dim, nexist);
		if (last[j] < i) {
		    Inner_Product(P->Constraint[j]+1, all_val->p, P->Dimension+1, 
				  &tmp);
		    if (value_neg_p(tmp))
			goto out2;
		    if (j < P->NbEq && value_pos_p(tmp))
			goto out2;
		}
	    }
	    --alternate;
	    --i;
	    continue;
	}
	assert(max_set && min_set);
	if (value_lt(max, min))
	    goto out2;
	if (value_ne(max, min)) {
	    Matrix *M = Matrix_Alloc(dim+i, 1+P->Dimension+1);
	    for (int j = 0; j < dim+i; ++j) {
		value_set_si(M->p[j][1+j], -1);
		value_assign(M->p[j][1+P->Dimension], all_val->p[j]);
	    }
	    Polyhedron *Q = AddConstraints(M->p[0], dim+i, P, MaxRays);
	    Matrix_Free(M);
	    Q = DomainConstraintSimplify(Q, MaxRays);
	    Vector *sample = Polyhedron_Sample(Q, MaxRays);
	    in = !!sample;
	    if (sample)
		Vector_Free(sample);
	    Polyhedron_Free(Q);
	    goto out2;
	}
	assert(value_eq(max, min));
	value_assign(all_val->p[dim+i], max);
	alternate = nexist - 1;
    }
    in = true;
out2:
    value_clear(min);
    value_clear(max);
out:
    Vector_Free(all_val);
    value_clear(tmp);
    return in || (P->next && in_domain(P->next, val, dim, MaxRays));
}

void compute_evalue(evalue *e, Value *val, Value *res)
{
    double d = compute_evalue(e, val);
    if (d > 0)
	d += .25;
    else
	d -= .25;
    value_set_double(*res, d);
}

Vector *max_term::eval(Value *val, unsigned MaxRays) const
{
    if (dim == domain->D->Dimension) {
	if (!in_domain(domain->D, val))
	    return NULL;
    } else {
	if (!in_domain(domain->D, val, dim, MaxRays))
	    return NULL;
    }
    Vector *res = Vector_Alloc(max.size());
    for (int i = 0; i < max.size(); ++i) {
	compute_evalue(max[i], val, &res->p[i]);
    }
    return res;
}

static Matrix *remove_equalities(Polyhedron **P, unsigned nparam, unsigned MaxRays);

Vector *Polyhedron_not_empty(Polyhedron *P, unsigned MaxRays)
{
    Polyhedron *Porig = P;
    Vector *sample = NULL;

    POL_ENSURE_VERTICES(P);
    if (emptyQ2(P))
	return NULL;

    for (int i = 0; i < P->NbRays; ++i)
	if (value_one_p(P->Ray[i][1+P->Dimension])) {
	    sample = Vector_Alloc(P->Dimension + 1);
	    Vector_Copy(P->Ray[i]+1, sample->p, P->Dimension+1);
	    return sample;
	}

    Matrix *T = remove_equalities(&P, 0, MaxRays);
    if (P)
	sample = Polyhedron_Sample(P, MaxRays);
    if (sample) {
	if (T) {
	    Vector *P_sample = Vector_Alloc(Porig->Dimension + 1);
	    Matrix_Vector_Product(T, sample->p, P_sample->p);
	    Vector_Free(sample);
	    sample = P_sample;
	}
    }
    if (T) {
	Polyhedron_Free(P);
	Matrix_Free(T);
    }

    return sample;
}

struct split {
    evalue *constraint;
    enum sign { le, ge, lge } sign;

    split (evalue *c, enum sign s) : constraint(c), sign(s) {}
};

static void split_on(const split& sp, EDomain *D, 
		     EDomain **Dlt, EDomain **Deq, EDomain **Dgt,
		     barvinok_options *options)
{
    Matrix *M, *M2;
    EDomain *ED[3];
    Polyhedron *D2;
    Value mone;
    value_init(mone);
    value_set_si(mone, -1);
    *Dlt = NULL;
    *Deq = NULL;
    *Dgt = NULL;
    vector<EDomain_floor *> new_floors;
    M = D->add_ge_constraint(sp.constraint, new_floors);
    if (sp.sign == split::lge || sp.sign == split::ge) {
	M2 = Matrix_Copy(M);
	value_decrement(M2->p[M2->NbRows-1][M2->NbColumns-1],
			M2->p[M2->NbRows-1][M2->NbColumns-1]);
	D2 = Constraints2Polyhedron(M2, options->MaxRays);
	ED[2] = new EDomain(D2, D, new_floors);
	Polyhedron_Free(D2);
	Matrix_Free(M2);
    } else
	ED[2] = NULL;
    if (sp.sign == split::lge || sp.sign == split::le) {
	M2 = Matrix_Copy(M);
	Vector_Scale(M2->p[M2->NbRows-1]+1, M2->p[M2->NbRows-1]+1,
		     mone, M2->NbColumns-1);
	value_decrement(M2->p[M2->NbRows-1][M2->NbColumns-1],
			M2->p[M2->NbRows-1][M2->NbColumns-1]);
	D2 = Constraints2Polyhedron(M2, options->MaxRays);
	ED[0] = new EDomain(D2, D, new_floors);
	Polyhedron_Free(D2);
	Matrix_Free(M2);
    } else
	ED[0] = NULL;

    assert(sp.sign == split::lge || sp.sign == split::ge || sp.sign == split::le);
    value_set_si(M->p[M->NbRows-1][0], 0);
    D2 = Constraints2Polyhedron(M, options->MaxRays);
    ED[1] = new EDomain(D2, D, new_floors);
    Polyhedron_Free(D2);
    Matrix_Free(M);

    Vector *sample = D->sample;
    if (sample && new_floors.size() > 0) {
	assert(sample->Size == D->D->Dimension+1);
	sample = Vector_Alloc(D->D->Dimension+new_floors.size()+1);
	Vector_Copy(D->sample->p, sample->p, D->D->Dimension);
	value_set_si(sample->p[D->D->Dimension+new_floors.size()], 1);
	for (int i = 0; i < new_floors.size(); ++i)
	    compute_evalue(new_floors[i]->e, sample->p, sample->p+D->D->Dimension+i);
    }

    for (int i = 0; i < new_floors.size(); ++i)
	EDomain_floor::unref(new_floors[i]);

    for (int i = 0; i < 3; ++i) {
	if (!ED[i])
	    continue;
	if (sample &&
		in_domain(ED[i]->D, sample->p, sample->Size-1, options->MaxRays)) {
	    ED[i]->sample = Vector_Alloc(sample->Size);
	    Vector_Copy(sample->p, ED[i]->sample->p, sample->Size);
	} else if (emptyQ2(ED[i]->D) ||
		    (options->emptiness_check == 1 &&
		     !(ED[i]->sample = Polyhedron_not_empty(ED[i]->D,
							    options->MaxRays)))) {
	    delete ED[i];
	    ED[i] = NULL;
	}
    }
    *Dlt = ED[0];
    *Deq = ED[1];
    *Dgt = ED[2];
    value_clear(mone);
    if (sample != D->sample)
	Vector_Free(sample);
}

ostream & operator<< (ostream & os, const vector<int> & v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        if (i)
            os << ", ";
        os << v[i]; 
    }
    os << "]";
    return os;
}

static bool isTranslation(Matrix *M)
{
    unsigned i, j;
    if (M->NbRows != M->NbColumns)
	return False;

    for (i = 0;i < M->NbRows; i ++)
	for (j = 0; j < M->NbColumns-1; j ++)
	    if (i == j) {
		if(value_notone_p(M->p[i][j]))
		    return False;
	    } else {
		if(value_notzero_p(M->p[i][j]))
		    return False;
	    }
    return value_one_p(M->p[M->NbRows-1][M->NbColumns-1]);
}

static Matrix *compress_parameters(Polyhedron **P, Polyhedron **C,
				   unsigned nparam, unsigned MaxRays)
{
    Matrix *M, *T, *CP;

    /* compress_parms doesn't like equalities that only involve parameters */
    for (int i = 0; i < (*P)->NbEq; ++i)
	assert(First_Non_Zero((*P)->Constraint[i]+1, (*P)->Dimension-nparam) != -1);

    M = Matrix_Alloc((*P)->NbEq, (*P)->Dimension+2);
    Vector_Copy((*P)->Constraint[0], M->p[0], (*P)->NbEq * ((*P)->Dimension+2));
    CP = compress_parms(M, nparam);
    Matrix_Free(M);

    if (isTranslation(CP)) {
	Matrix_Free(CP);
	return NULL;
    }

    T = align_matrix(CP, (*P)->Dimension+1);
    *P = Polyhedron_Preimage(*P, T, MaxRays);
    Matrix_Free(T);

    *C = Polyhedron_Preimage(*C, CP, MaxRays);

    return CP;
}

static Matrix *remove_equalities(Polyhedron **P, unsigned nparam, unsigned MaxRays)
{
    /* Matrix "view" of equalities */
    Matrix M;
    M.NbRows = (*P)->NbEq;
    M.NbColumns = (*P)->Dimension+2;
    M.p_Init = (*P)->p_Init;
    M.p = (*P)->Constraint;

    Matrix *T = compress_variables(&M, nparam);

    if (!T) {
	*P = NULL;
	return NULL;
    }
    if (isIdentity(T)) {
	Matrix_Free(T);
	T = NULL;
    } else
	*P = Polyhedron_Preimage(*P, T, MaxRays);

    return T;
}

void construct_rational_vertices(Param_Polyhedron *PP, Matrix *T, unsigned dim, 
				 int nparam, vector<indicator_term *>& vertices)
{
    int i;
    Param_Vertices *PV;
    Value lcm, tmp;
    value_init(lcm);
    value_init(tmp);

    vec_ZZ v;
    v.SetLength(nparam+1);

    evalue factor;
    value_init(factor.d);
    value_init(factor.x.n);
    value_set_si(factor.x.n, 1);
    value_set_si(factor.d, 1);

    for (i = 0, PV = PP->V; PV; ++i, PV = PV->next) {
	indicator_term *term = new indicator_term(dim, i);
	vertices.push_back(term);
	Matrix *M = Matrix_Alloc(PV->Vertex->NbRows+nparam+1, nparam+1);
	value_set_si(lcm, 1);
	for (int j = 0; j < PV->Vertex->NbRows; ++j)
	    value_lcm(lcm, PV->Vertex->p[j][nparam+1], &lcm);
	value_assign(M->p[M->NbRows-1][M->NbColumns-1], lcm);
	for (int j = 0; j < PV->Vertex->NbRows; ++j) {
	    value_division(tmp, lcm, PV->Vertex->p[j][nparam+1]);
	    Vector_Scale(PV->Vertex->p[j], M->p[j], tmp, nparam+1);
	}
	for (int j = 0; j < nparam; ++j)
	    value_assign(M->p[PV->Vertex->NbRows+j][j], lcm);
	if (T) {
	    Matrix *M2 = Matrix_Alloc(T->NbRows, M->NbColumns);
	    Matrix_Product(T, M, M2);
	    Matrix_Free(M);
	    M = M2;
	}
	for (int j = 0; j < dim; ++j) {
	    values2zz(M->p[j], v, nparam+1);
	    term->vertex[j] = multi_monom(v);
	    value_assign(factor.d, lcm);
	    emul(&factor, term->vertex[j]);
	}
	Matrix_Free(M);
    }
    assert(i == PP->nbV);
    free_evalue_refs(&factor);
    value_clear(lcm);
    value_clear(tmp);
}

/* An auxiliary class that keeps a reference to an evalue
 * and frees it when it goes out of scope.
 */
struct temp_evalue {
    evalue *E;
    temp_evalue() : E(NULL) {}
    temp_evalue(evalue *e) : E(e) {}
    operator evalue* () const { return E; }
    evalue *operator=(evalue *e) {
	if (E) {
	    free_evalue_refs(E); 
	    delete E;
	}
	E = e;
	return E;
    }
    ~temp_evalue() {
	if (E) {
	    free_evalue_refs(E); 
	    delete E;
	}
    }
};

static vector<max_term*> lexmin(indicator& ind, unsigned nparam,
				    vector<int> loc)
{
    vector<max_term*> maxima;
    map<indicator_term *, int >::iterator i;
    vector<int> best_score;
    vector<int> second_score;
    vector<int> neg_score;

    do {
	indicator_term *best = NULL, *second = NULL, *neg = NULL,
		       *neg_eq = NULL, *neg_le = NULL;
	for (i = ind.order.pred.begin(); i != ind.order.pred.end(); ++i) {
	    vector<int> score;
	    if ((*i).second != 0)
		continue;
	    indicator_term *term = (*i).first;
	    if (term->sign == 0) {
		ind.expand_rational_vertex(term);
		break;
	    }

	    if (ind.order.eq.find(term) != ind.order.eq.end()) {
		int j;
		if (ind.order.eq[term].size() <= 1)
		    continue;
		for (j = 1; j < ind.order.eq[term].size(); ++j)
		    if (ind.order.pred[ind.order.eq[term][j]] != 0)
			break;
		if (j < ind.order.eq[term].size())
		    continue;
		score.push_back(ind.order.eq[term].size());
	    } else
		score.push_back(0);
	    if (ind.order.le.find(term) != ind.order.le.end())
		score.push_back(ind.order.le[term].size());
	    else
		score.push_back(0);
	    if (ind.order.lt.find(term) != ind.order.lt.end())
		score.push_back(-ind.order.lt[term].size());
	    else
		score.push_back(0);

	    if (term->sign > 0) {
		if (!best || score < best_score) {
		    second = best;
		    second_score = best_score;
		    best = term;
		    best_score = score;
		} else if (!second || score < second_score) {
		    second = term;
		    second_score = score;
		}
	    } else {
		if (!neg_eq && ind.order.eq.find(term) != ind.order.eq.end()) {
		    for (int j = 1; j < ind.order.eq[term].size(); ++j)
			if (ind.order.eq[term][j]->sign != term->sign) {
			    neg_eq = term;
			    break;
			}
		}
		if (!neg_le && ind.order.le.find(term) != ind.order.le.end())
		    neg_le = term;
		if (!neg || score < neg_score) {
		    neg = term;
		    neg_score = score;
		}
	    }
	}
	if (i != ind.order.pred.end())
	    continue;

	if (!best && neg_eq) {
	    assert(ind.order.eq[neg_eq].size() != 0);
	    bool handled = ind.handle_equal_numerators(neg_eq);
	    assert(handled);
	    continue;
	}

	if (!best && neg_le) {
	    /* The smallest term is negative and <= some positive term */
	    best = neg_le;
	    neg = NULL;
	}

	if (!best) {
	    /* apparently there can be negative initial term on empty domains */
	    if (ind.options->emptiness_check == 1)
		assert(!neg);
	    break;
	}

	if (!second && !neg) {
	    indicator_term *rat = NULL;
	    assert(best);
	    if (ind.order.le[best].size() == 0) {
		if (ind.order.eq[best].size() != 0) {
		    bool handled = ind.handle_equal_numerators(best);
		    if (ind.options->emptiness_check == 1)
			assert(handled);
		    /* If !handled then the leading coefficient is bigger than one;
		     * must be an empty domain
		     */
		    if (handled)
			continue;
		    else
			break;
		}
		maxima.push_back(ind.create_max_term(best));
		break;
	    }
	    for (int j = 0; j < ind.order.le[best].size(); ++j) {
		if (ind.order.le[best][j]->sign == 0) {
		    if (!rat && ind.order.pred[ind.order.le[best][j]] == 1)
			rat = ind.order.le[best][j];
		} else if (ind.order.le[best][j]->sign > 0) {
		    second = ind.order.le[best][j];
		    break;
		} else if (!neg)
		    neg = ind.order.le[best][j];
	    }

	    if (!second && !neg) {
		assert(rat);
		ind.order.unset_le(best, rat);
		ind.expand_rational_vertex(rat);
		continue;
	    }

	    if (!second)
		second = neg;

	    ind.order.unset_le(best, second);
	}

	if (!second)
	    second = neg;

	unsigned dim = best->den.NumCols();
	temp_evalue diff;
	order_sign sign;
	for (int k = 0; k < dim; ++k) {
	    diff = ediff(best->vertex[k], second->vertex[k]);
	    sign = evalue_sign(diff, ind.D, ind.options->MaxRays);

	    /* neg can never be smaller than best, unless it may still cancel
	     * This optimization is disabled for now since we can't know
	     * for sure that neg won't cancel.  It may not have been compared
	     * yet to another negative term that is equal of less than or eqaul
	     * to some positive terms.
	     */
	    if (0 && second == neg &&
		ind.order.eq.find(neg) == ind.order.eq.end() &&
		ind.order.le.find(neg) == ind.order.le.end()) {
		if (sign == order_ge)
		    sign = order_eq;
		if (sign == order_unknown)
		    sign = order_le;
	    }

	    if (sign != order_eq)
		break;
	    if (!EVALUE_IS_ZERO(*diff))
		ind.substitute(diff);
	}
	if (sign == order_eq) {
	    ind.order.set_equal(best, second);
	    continue;
	}
	if (sign == order_lt) {
	    ind.order.lt[best].push_back(second);
	    ind.order.pred[second]++;
	    continue;
	}
	if (sign == order_gt) {
	    ind.order.lt[second].push_back(best);
	    ind.order.pred[best]++;
	    continue;
	}

	split sp(diff, sign == order_le ? split::le :
		       sign == order_ge ? split::ge : split::lge);

	EDomain *Dlt, *Deq, *Dgt;
	split_on(sp, ind.D, &Dlt, &Deq, &Dgt, ind.options);
	if (ind.options->emptiness_check == 1)
	    assert(Dlt || Deq || Dgt);
	else if (!(Dlt || Deq || Dgt))
	    /* Must have been empty all along */
	    break;

	if (Deq && (Dlt || Dgt)) {
	    int locsize = loc.size();
	    loc.push_back(0);
	    indicator indeq(ind, Deq);
	    Deq = NULL;
	    indeq.substitute(diff);
	    vector<max_term*> maxeq = lexmin(indeq, nparam, loc);
	    maxima.insert(maxima.end(), maxeq.begin(), maxeq.end());
	    loc.resize(locsize);
	}
	if (Dgt && Dlt) {
	    int locsize = loc.size();
	    loc.push_back(1);
	    indicator indgt(ind, Dgt);
	    Dgt = NULL;
	    /* we don't know the new location of these terms in indgt */
	    /*
	    indgt.order.lt[second].push_back(best);
	    indgt.order.pred[best]++;
	    */
	    vector<max_term*> maxgt = lexmin(indgt, nparam, loc);
	    maxima.insert(maxima.end(), maxgt.begin(), maxgt.end());
	    loc.resize(locsize);
	}

	if (Deq) {
	    loc.push_back(0);
	    ind.substitute(diff);
	    ind.set_domain(Deq);
	}
	if (Dlt) {
	    loc.push_back(-1);
	    ind.order.lt[best].push_back(second);
	    ind.order.pred[second]++;
	    ind.set_domain(Dlt);
	}
	if (Dgt) {
	    loc.push_back(1);
	    ind.order.lt[second].push_back(best);
	    ind.order.pred[best]++;
	    ind.set_domain(Dgt);
	}
    } while(1);

    return maxima;
}

static vector<max_term*> lexmin(Polyhedron *P, Polyhedron *C,
				barvinok_options *options)
{
    unsigned nparam = C->Dimension;
    Param_Polyhedron *PP = NULL;
    Polyhedron *CEq = NULL, *pVD;
    Matrix *CT = NULL;
    Matrix *T = NULL, *CP = NULL;
    Param_Domain *D, *next;
    Param_Vertices *V;
    Polyhedron *Porig = P;
    Polyhedron *Corig = C;
    vector<max_term*> all_max;
    Polyhedron *Q;
    unsigned P2PSD_MaxRays;

    if (emptyQ2(P))
	return all_max;

    POL_ENSURE_VERTICES(P);

    if (emptyQ2(P))
	return all_max;

    assert(P->NbBid == 0);

    if (P->NbEq > 0) {
	if (nparam > 0)
	    CP = compress_parameters(&P, &C, nparam, options->MaxRays);
	Q = P;
	T = remove_equalities(&P, nparam, options->MaxRays);
	if (P != Q && Q != Porig)
	    Polyhedron_Free(Q);
	if (!P) {
	    if (C != Corig)
		Polyhedron_Free(C);
	    return all_max;
	}
    }

    if (options->MaxRays & POL_NO_DUAL)
	P2PSD_MaxRays = 0;
    else
	P2PSD_MaxRays = options->MaxRays;

    Q = P;
    PP = Polyhedron2Param_SimplifiedDomain(&P, C, P2PSD_MaxRays, &CEq, &CT);
    if (P != Q && Q != Porig)
	Polyhedron_Free(Q);

    if (CT) {
	if (isIdentity(CT)) {
	    Matrix_Free(CT);
	    CT = NULL;
	} else
	    nparam = CT->NbRows - 1;
    }
    assert(!CT);

    unsigned dim = P->Dimension - nparam;

    int nd;
    for (nd = 0, D=PP->D; D; ++nd, D=D->next);
    Polyhedron **fVD = new Polyhedron*[nd];

    indicator_constructor ic(P, dim, PP, T);

    vector<indicator_term *> all_vertices;
    construct_rational_vertices(PP, T, T ? T->NbRows-nparam-1 : dim,
				nparam, all_vertices);

    for (nd = 0, D=PP->D; D; D=next) {
	next = D->next;

	Polyhedron *rVD = reduce_domain(D->Domain, CT, CEq,
					fVD, nd, options->MaxRays);
	if (!rVD)
	    continue;

	pVD = CT ? DomainImage(rVD,CT,options->MaxRays) : rVD;

	EDomain *epVD = new EDomain(pVD);
	indicator ind(ic, D, epVD, options);

	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    ind.add(all_vertices[_i]);
	END_FORALL_PVertex_in_ParamPolyhedron;

	ind.remove_initial_rational_vertices();

	vector<int> loc;
	vector<max_term*> maxima = lexmin(ind, nparam, loc);
	if (CP)
	    for (int j = 0; j < maxima.size(); ++j)
		maxima[j]->substitute(CP, options->MaxRays);
	all_max.insert(all_max.end(), maxima.begin(), maxima.end());

	++nd;
	if (rVD != pVD)
	    Domain_Free(pVD);
	Domain_Free(rVD);
    }
    for (int i = 0; i < all_vertices.size(); ++i)
	delete all_vertices[i];
    if (CP)
	Matrix_Free(CP);
    if (T)
	Matrix_Free(T);
    Param_Polyhedron_Free(PP);
    if (CEq)
	Polyhedron_Free(CEq);
    for (--nd; nd >= 0; --nd) {
	Domain_Free(fVD[nd]);
    }
    delete [] fVD;
    if (C != Corig)
	Polyhedron_Free(C);
    if (P != Porig)
	Polyhedron_Free(P);

    return all_max;
}

static void verify_results(Polyhedron *A, Polyhedron *C, 
			   vector<max_term*>& maxima, int m, int M,
			   int print_all, unsigned MaxRays);

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *MA;
    int bignum;
    char **iter_names, **param_names;
    int c, ind = 0;
    int range = 0;
    int verify = 0;
    int print_all = 0;
    int m = INT_MAX, M = INT_MIN, r;
    int print_solution = 1;
    struct barvinok_options *options;

    options = barvinok_options_new_with_defaults();

    while ((c = getopt_long(argc, argv, "TAm:M:r:V", lexmin_options, &ind)) != -1) {
	switch (c) {
	case NO_EMPTINESS_CHECK:
	    options->emptiness_check = 0;
	    break;
	case 'T':
	    verify = 1;
	    break;
	case 'A':
	    print_all = 1;
	    break;
	case 'm':
	    m = atoi(optarg);
	    verify = 1;
	    break;
	case 'M':
	    M = atoi(optarg);
	    verify = 1;
	    break;
	case 'r':
	    M = atoi(optarg);
	    m = -M;
	    verify = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    MA = Matrix_Read();
    C = Constraints2Polyhedron(MA, options->MaxRays);
    Matrix_Free(MA);
    fscanf(stdin, " %d", &bignum);
    assert(bignum == -1);
    MA = Matrix_Read();
    A = Constraints2Polyhedron(MA, options->MaxRays);
    Matrix_Free(MA);

    if (A->Dimension >= VBIGDIM)
	r = VSRANGE;
    else if (A->Dimension >= BIGDIM)
	r = SRANGE;
    else
	r = RANGE;
    if (M == INT_MIN)
	M = r;
    if (m == INT_MAX)
	m = -r;

    if (verify && m > M) {
	fprintf(stderr,"Nothing to do: min > max !\n");
	return(0);
    }
    if (verify)
	print_solution = 0;

    iter_names = util_generate_names(A->Dimension - C->Dimension, "i");
    param_names = util_generate_names(C->Dimension, "p");
    if (print_solution) {
	Polyhedron_Print(stdout, P_VALUE_FMT, A);
	Polyhedron_Print(stdout, P_VALUE_FMT, C);
    }
    vector<max_term*> maxima = lexmin(A, C, options);
    if (print_solution)
	for (int i = 0; i < maxima.size(); ++i)
	    maxima[i]->print(cout, param_names);

    if (verify)
	verify_results(A, C, maxima, m, M, print_all, options->MaxRays);

    for (int i = 0; i < maxima.size(); ++i)
	delete maxima[i];

    util_free_names(A->Dimension - C->Dimension, iter_names);
    util_free_names(C->Dimension, param_names);
    Polyhedron_Free(A);
    Polyhedron_Free(C);

    free(options);

    return 0;
}

static bool lexmin(int pos, Polyhedron *P, Value *context)
{
    Value LB, UB, k;
    int lu_flags;
    bool found = false;

    if (emptyQ(P))
	return false;

    value_init(LB); value_init(UB); value_init(k);
    value_set_si(LB,0);
    value_set_si(UB,0);
    lu_flags = lower_upper_bounds(pos,P,context,&LB,&UB);
    assert(!(lu_flags & LB_INFINITY));

    value_set_si(context[pos],0);
    if (!lu_flags && value_lt(UB,LB)) {
        value_clear(LB); value_clear(UB); value_clear(k);
	return false;
    }
    if (!P->next) {
	value_assign(context[pos], LB);
        value_clear(LB); value_clear(UB); value_clear(k);
	return true;
    }
    for (value_assign(k,LB); lu_flags || value_le(k,UB); value_increment(k,k)) {
        value_assign(context[pos],k);
	if ((found = lexmin(pos+1, P->next, context)))
	    break;
    }
    if (!found)
	value_set_si(context[pos],0);
    value_clear(LB); value_clear(UB); value_clear(k);
    return found;
}

static void print_list(FILE *out, Value *z, char* brackets, int len)
{
    fprintf(out, "%c", brackets[0]);
    value_print(out, VALUE_FMT,z[0]);
    for (int k = 1; k < len; ++k) {
	fprintf(out, ", ");
	value_print(out, VALUE_FMT,z[k]);
    }
    fprintf(out, "%c", brackets[1]);
}

static int check_poly(Polyhedron *S, Polyhedron *CS, vector<max_term*>& maxima, 
		      int nparam, int pos, Value *z, int print_all, int st,
		      unsigned MaxRays)
{
    if (pos == nparam) {
	int k;
	bool found = lexmin(1, S, z);

	if (print_all) {
	    printf("lexmin");
	    print_list(stdout, z+S->Dimension-nparam+1, "()", nparam);
	    printf(" = ");
	    if (found)
		print_list(stdout, z+1, "[]", S->Dimension-nparam);
	    printf(" ");
	}

	Vector *min = NULL;
	for (int i = 0; i < maxima.size(); ++i)
	    if ((min = maxima[i]->eval(z+S->Dimension-nparam+1, MaxRays)))
		break;

	int ok = !(found ^ !!min);
	if (found && min)
	    for (int i = 0; i < S->Dimension-nparam; ++i)
		if (value_ne(z[1+i], min->p[i])) {
		    ok = 0;
		    break;
		}
	if (!ok) {
	    printf("\n"); 
	    fflush(stdout);
	    fprintf(stderr, "Error !\n");
	    fprintf(stderr, "lexmin");
	    print_list(stderr, z+S->Dimension-nparam+1, "()", nparam);
	    fprintf(stderr, " should be ");
	    if (found)
		print_list(stderr, z+1, "[]", S->Dimension-nparam);
	    fprintf(stderr, " while digging gives ");
	    if (min)
		print_list(stderr, min->p, "[]", S->Dimension-nparam);
	    fprintf(stderr, ".\n");
	    return 0;
	} else if (print_all)
	    printf("OK.\n");
	if (min)
	    Vector_Free(min);

	for (k = 1; k <= S->Dimension-nparam; ++k)
	    value_set_si(z[k], 0);
    } else {
	Value tmp;
	Value LB, UB;
	value_init(tmp);
	value_init(LB);
	value_init(UB);
	int ok = 
	    !(lower_upper_bounds(1+pos, CS, &z[S->Dimension-nparam], &LB, &UB));
	for (value_assign(tmp,LB); value_le(tmp,UB); value_increment(tmp,tmp)) {
	    if (!print_all) {
		int k = VALUE_TO_INT(tmp);
		if (!pos && !(k%st)) {
		    printf("o");
		    fflush(stdout);
		}
	    }
	    value_assign(z[pos+S->Dimension-nparam+1],tmp);
	    if (!check_poly(S, CS->next, maxima, nparam, pos+1, z, print_all, st,
			    MaxRays)) {
		value_clear(tmp);
		value_clear(LB);
		value_clear(UB);
		return 0;
	    }
	}
	value_set_si(z[pos+S->Dimension-nparam+1],0);
	value_clear(tmp);
	value_clear(LB);
	value_clear(UB);
    }
    return 1;
}

void verify_results(Polyhedron *A, Polyhedron *C, vector<max_term*>& maxima, 
		    int m, int M, int print_all, unsigned MaxRays)
{
    Polyhedron *CC, *CC2, *CS, *S;
    unsigned nparam = C->Dimension;
    Value *p;
    int i;
    int st;

    CC = Polyhedron_Project(A, nparam);
    CC2 = DomainIntersection(C, CC, MaxRays);
    Domain_Free(CC);
    CC = CC2;

    /* Intersect context with range */
    if (nparam > 0) {
	Matrix *MM;
	Polyhedron *U;

	MM = Matrix_Alloc(2*C->Dimension, C->Dimension+2);
	for (int i = 0; i < C->Dimension; ++i) {
	    value_set_si(MM->p[2*i][0], 1);
	    value_set_si(MM->p[2*i][1+i], 1);
	    value_set_si(MM->p[2*i][1+C->Dimension], -m);
	    value_set_si(MM->p[2*i+1][0], 1);
	    value_set_si(MM->p[2*i+1][1+i], -1);
	    value_set_si(MM->p[2*i+1][1+C->Dimension], M);
	}
	CC2 = AddConstraints(MM->p[0], 2*CC->Dimension, CC, MaxRays);
	U = Universe_Polyhedron(0);
	CS = Polyhedron_Scan(CC2, U, MaxRays & POL_NO_DUAL ? 0 : MaxRays);
	Polyhedron_Free(U);
	Polyhedron_Free(CC2);
	Matrix_Free(MM);
    } else
	CS = NULL;

    p = ALLOCN(Value, A->Dimension+2);
    for (i=0; i <= A->Dimension; i++) {
	value_init(p[i]);
	value_set_si(p[i],0);
    }
    value_init(p[i]);
    value_set_si(p[i], 1);

    S = Polyhedron_Scan(A, C, MaxRays & POL_NO_DUAL ? 0 : MaxRays);

    if (!print_all && C->Dimension > 0) {
	if (M-m > 80)
	    st = 1+(M-m)/80;
	else
	    st = 1;
	for (int i = m; i <= M; i += st)
	    printf(".");
	printf( "\r" );
	fflush(stdout);
    }

    if (S) {
	if (!(CS && emptyQ2(CS)))
	    check_poly(S, CS, maxima, nparam, 0, p, print_all, st, MaxRays);
	Domain_Free(S);
    }

    if (!print_all)
	printf("\n");

    for (i=0; i <= (A->Dimension+1); i++) 
	value_clear(p[i]);
    free(p);
    if (CS)
	Domain_Free(CS);
    Polyhedron_Free(CC);
}
