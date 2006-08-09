#include <iostream>
#include <vector>
#include <gmp.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/barvinok.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "decomposer.h"
#include "lattice_point.h"
#include "reduce_domain.h"
#include "mat_util.h"
#include "combine.h"
#include "sample.h"
#include "fdstream.h"
#include "config.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    (POL_NO_DUAL | POL_INTEGER)
#else
#define MAXRAYS  600
#endif

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
struct option options[] = {
    { "verify",     no_argument,  0,  'T' },
    { "print-all",  no_argument,  0,  'A' },
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

static void evalue_denom(evalue *e, Value *d)
{
    if (value_notzero_p(e->d)) {
	value_lcm(*d, e->d, d);
	return;
    }
    int offset = type_offset(e->x.p);
    for (int i = e->x.p->size-1; i >= offset; --i)
	evalue_denom(&e->x.p->arr[i], d);
}

static void evalue_print(std::ostream& o, evalue *e, char **p);
static void evalue_print(std::ostream& o, evalue *e, char **p, int d)
{
    if (value_notzero_p(e->d)) {
	o << VALUE_TO_INT(e->x.n) * (d / VALUE_TO_INT(e->d));
	return;
    }
    assert(e->x.p->type == polynomial || e->x.p->type == flooring ||
	   e->x.p->type == fractional);
    int offset = type_offset(e->x.p);
    for (int i = e->x.p->size-1; i >= offset; --i) {
	if (EVALUE_IS_ZERO(e->x.p->arr[i]))
	    continue;
	if (i != e->x.p->size-1 && 
		(value_zero_p(e->x.p->arr[i].d) ||
		value_pos_p(e->x.p->arr[i].x.n)))
	    o << "+";
	if (i == offset || !(value_one_p(e->x.p->arr[i].x.n) && 
			     d == VALUE_TO_INT(e->x.p->arr[i].d))) {
	    if (value_zero_p(e->x.p->arr[i].d))
		o << "(";
	    evalue_print(o, &e->x.p->arr[i], p, d);
	    if (value_zero_p(e->x.p->arr[i].d))
		o << ")";
	    if (i != offset)
		o << "*";
	}
	for (int j = 0; j < i-offset; ++j) {
	    if (j != 0)
		o << "*";
	    if (e->x.p->type == flooring) {
		o << "[";
		evalue_print(o, &e->x.p->arr[0], p);
		o << "]";
	    } else if (e->x.p->type == fractional) {
		o << "{";
		evalue_print(o, &e->x.p->arr[0], p);
		o << "}";
	    } else
		o << p[e->x.p->pos-1];
	}
    }
}

static void evalue_print(std::ostream& o, evalue *e, char **p)
{
    Value d;
    value_init(d);
    value_set_si(d, 1);
    evalue_denom(e, &d);
    if (value_notone_p(d))
	o << "(";
    evalue_print(o, e, p, VALUE_TO_INT(d));
    if (value_notone_p(d))
	o << ")/" << VALUE_TO_INT(d);
    value_clear(d);
}

struct indicator_term {
    int sign;
    mat_ZZ den;
    evalue **vertex;

    indicator_term(unsigned dim) {
	den.SetDims(dim, dim);
	vertex = new evalue* [dim];
    }
    indicator_term(const indicator_term& src) {
	sign = src.sign;
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
    void print(ostream& os, char **p);
    void substitute(Matrix *T);
    void substitute(evalue *fract, evalue *val);
    void substitute(int pos, evalue *val);
    void reduce_in_domain(Polyhedron *D);
};

void indicator_term::reduce_in_domain(Polyhedron *D)
{
    for (int k = 0; k < den.NumCols(); ++k) {
	reduce_evalue_in_domain(vertex[k], D);
	if (evalue_range_reduction_in_domain(vertex[k], D))
	    reduce_evalue(vertex[k]);
    }
}

void indicator_term::print(ostream& os, char **p)
{
    unsigned dim = den.NumCols();
    unsigned factors = den.NumRows();
    if (sign > 0)
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
}

/* Perform the substitution specified by T on the variables.
 * T has dimension (newdim+1) x (olddim + nparam + 1).
 * The substitution is performed as in gen_fun::substitute
 */
void indicator_term::substitute(Matrix *T)
{
    unsigned dim = den.NumCols();
    unsigned newdim = T->NbRows - 1;
    unsigned nparam = T->NbColumns - dim - 1;
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

    indicator_constructor(Polyhedron *P, unsigned dim, unsigned nbV) :
		vertex_decomposer(P, nbV, this) {
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
    lex_order_rows(term->den);

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
}

void indicator_constructor::substitute(Matrix *T)
{
    for (int i = 0; i < nbV; ++i)
	for (int j = 0; j < terms[i].size(); ++j)
	    terms[i][j]->substitute(T);
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

struct indicator {
    vector<indicator_term*> term;

    indicator() {}
    indicator(const indicator& ind) {
	for (int i = 0; i < ind.term.size(); ++i)
	    term.push_back(new indicator_term(*ind.term[i]));
    }
    ~indicator() {
	for (int i = 0; i < term.size(); ++i)
	    delete term[i];
    }

    void print(ostream& os, char **p);
    void simplify();
    void peel(int i, int j);
    void combine(int i, int j);
    void substitute(evalue *equation);
    void reduce_in_domain(Polyhedron *D);
};

void indicator::reduce_in_domain(Polyhedron *D)
{
    for (int i = 0; i < term.size(); ++i)
	term[i]->reduce_in_domain(D);
}

void indicator::print(ostream& os, char **p)
{
    for (int i = 0; i < term.size(); ++i) {
	term[i]->print(os, p);
	os << endl;
    }
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

void indicator::combine(int i, int j)
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

    assert(n_i < 30);
    assert(n_j < 30);

    for (k = 0; k < (1 << n_i); ++k) {
	indicator_term *it = new indicator_term(*term[j]);
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
	term.push_back(it);
    }

    for (k = 0; k < (1 << n_j); ++k) {
	indicator_term *it = new indicator_term(*term[i]);
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
	term.push_back(it);
    }
    delete term[i];
    delete term[j];
    term.erase(term.begin()+j);
    term.erase(term.begin()+i);
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

static void add_coeff(Value *cons, int len, evalue *coeff, int pos)
{
    Value tmp;

    assert(value_notzero_p(coeff->d));

    value_init(tmp);

    value_lcm(cons[0], coeff->d, &tmp);
    value_division(tmp, tmp, cons[0]);
    Vector_Scale(cons, cons, tmp, len);
    value_division(tmp, cons[0], coeff->d);
    value_addmul(cons[pos], tmp, coeff->x.n);

    value_clear(tmp);
}

struct EDomain {
    Polyhedron		*D;
    Vector		*sample;
    vector<evalue *>	floors;

    EDomain(Polyhedron *D) {
	this->D = Polyhedron_Copy(D);
	sample = NULL;
    }
    EDomain(Polyhedron *D, vector<evalue *>floors) {
	this->D = Polyhedron_Copy(D);
	add_floors(floors);
	sample = NULL;
    }
    EDomain(Polyhedron *D, EDomain *ED, vector<evalue *>floors) {
	this->D = Polyhedron_Copy(D);
	add_floors(ED->floors);
	add_floors(floors);
	sample = NULL;
    }
    void add_floors(vector<evalue *>floors) {
	for (int i = 0; i < floors.size(); ++i) {
	    evalue *f = new evalue;
	    value_init(f->d);
	    evalue_copy(f, floors[i]);
	    this->floors.push_back(f);
	}
    }
    int find_floor(evalue *needle) {
	for (int i = 0; i < floors.size(); ++i)
	    if (eequal(needle, floors[i]))
		return i;
	return -1;
    }
    void print(FILE *out, char **p);
    ~EDomain() {
	for (int i = 0; i < floors.size(); ++i) {
	    free_evalue_refs(floors[i]);
	    delete floors[i];
	}
	Polyhedron_Free(D);
	if (sample)
	    Vector_Free(sample);
    }
};

void EDomain::print(FILE *out, char **p)
{
    fdostream os(dup(fileno(out)));
    for (int i = 0; i < floors.size(); ++i) {
	os << "floor " << i << ": [";
	evalue_print(os, floors[i], p);
	os << "]" << endl;
    }
    Polyhedron_Print(out, P_VALUE_FMT, D);
}

static int evalue2constraint_r(EDomain *D, evalue *E, Value *cons, int len);

static void add_fract(evalue *e, Value *cons, int len, int pos)
{
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    /* contribution of alpha * fract(X) is 
     *      alpha * X ...
     */
    assert(e->x.p->size == 3);
    evalue argument;
    value_init(argument.d);
    evalue_copy(&argument, &e->x.p->arr[0]);
    emul(&e->x.p->arr[2], &argument);
    evalue2constraint_r(NULL, &argument, cons, len);
    free_evalue_refs(&argument);

    /*	    - alpha * floor(X) */
    emul(&mone, &e->x.p->arr[2]);
    add_coeff(cons, len, &e->x.p->arr[2], pos);
    emul(&mone, &e->x.p->arr[2]);

    free_evalue_refs(&mone); 
}

static int evalue2constraint_r(EDomain *D, evalue *E, Value *cons, int len)
{
    int r = 0;
    if (value_zero_p(E->d) && E->x.p->type == fractional) {
	int i;
	assert(E->x.p->size == 3);
	r = evalue2constraint_r(D, &E->x.p->arr[1], cons, len);
	assert(value_notzero_p(E->x.p->arr[2].d));
	if (D && (i = D->find_floor(&E->x.p->arr[0])) >= 0) {
	    add_fract(E, cons, len, 1+D->D->Dimension-D->floors.size()+i);
	} else {
	    if (value_pos_p(E->x.p->arr[2].x.n)) {
		evalue coeff;
		value_init(coeff.d);
		value_init(coeff.x.n);
		value_set_si(coeff.d, 1);
		evalue_denom(&E->x.p->arr[0], &coeff.d);
		value_decrement(coeff.x.n, coeff.d);
		emul(&E->x.p->arr[2], &coeff);
		add_coeff(cons, len, &coeff, len-1);
		free_evalue_refs(&coeff);
	    }
	    r = 1;
	}
    } else if (value_zero_p(E->d)) {
	assert(E->x.p->type == polynomial);
	assert(E->x.p->size == 2);
	r = evalue2constraint_r(D, &E->x.p->arr[0], cons, len) || r;
	add_coeff(cons, len, &E->x.p->arr[1], E->x.p->pos);
    } else {
	add_coeff(cons, len, E, len-1);
    }
    return r;
}

static int evalue2constraint(EDomain *D, evalue *E, Value *cons, int len)
{
    Vector_Set(cons, 0, len);
    value_set_si(cons[0], 1);
    return evalue2constraint_r(D, E, cons, len);
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
    Polyhedron *domain;
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
	Polyhedron_Free(domain);
    }
};

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
    ::resolve_existential_vars(domain, dim);
}

void max_term::print(ostream& os, char **p) const
{
    char **names = p;
    if (dim < domain->Dimension) {
	resolve_existential_vars();
	names = new char * [domain->Dimension];
	int i;
	for (i = 0; i < dim; ++i)
	    names[i] = p[i];
	for ( ; i < domain->Dimension; ++i) {
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
    if (dim < domain->Dimension) {
	os << "exists ";
	print_varlist(os, domain->Dimension-dim, names+dim);
	os << " : ";
    }
    for (int i = 0; i < domain->NbConstraints; ++i) {
	int first = 1;
	int v = First_Non_Zero(domain->Constraint[i]+1, domain->Dimension);
	if (v == -1)
	    continue;
	if (i)
	    os << " && ";
	if (value_pos_p(domain->Constraint[i][v+1])) {
	    print_term(os, domain->Constraint[i][v+1], v, domain->Dimension,
		       names, NULL);
	    if (value_zero_p(domain->Constraint[i][0]))
		os << " = ";
	    else
		os << " >= ";
	    for (int j = v+1; j <= domain->Dimension; ++j) {
		value_oppose(tmp, domain->Constraint[i][1+j]);
		print_term(os, tmp, j, domain->Dimension,
			   names, &first);
	    }
	} else {
	    value_oppose(tmp, domain->Constraint[i][1+v]);
	    print_term(os, tmp, v, domain->Dimension,
		       names, NULL);
	    if (value_zero_p(domain->Constraint[i][0]))
		os << " = ";
	    else
		os << " <= ";
	    for (int j = v+1; j <= domain->Dimension; ++j)
		print_term(os, domain->Constraint[i][1+j], j, domain->Dimension,
			   names, &first);
	}
    }
    os << " }" << endl;
    value_clear(tmp);

    if (dim < domain->Dimension) {
	for (int i = dim; i < domain->Dimension; ++i)
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

/* T maps the compressed parameters to the original parameters,
 * while this max_term is based on the compressed parameters
 * and we want get the original parameters back.
 */
void max_term::substitute(Matrix *T, unsigned MaxRays)
{
    int nexist = 0;
    for (int i = 0; i < T->NbRows-1; ++i)
	if (value_notone_p(T->p[i][i]))
	    ++nexist;

    Matrix *M = Matrix_Alloc(T->NbRows + nexist, T->NbColumns);
    nexist = 0;
    for (int i = 0; i < T->NbRows-1; ++i) {
	Vector_Copy(T->p[i], M->p[i], T->NbColumns);
	if (value_notone_p(T->p[i][i]))
	    value_set_si(M->p[T->NbRows-1 + nexist++][i], 1);
    }
    value_assign(M->p[M->NbRows-1][M->NbColumns-1],
		 T->p[T->NbRows-1][T->NbColumns-1]);

    Polyhedron *D = DomainImage(domain, M, MaxRays);
    Polyhedron_Free(domain);
    domain = D;
    Matrix_Free(M);

    assert(T->NbRows == T->NbColumns);
    Matrix *inv = Matrix_Alloc(T->NbColumns, T->NbRows);
    int ok = Matrix_Inverse(T, inv);
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
    if (dim == domain->Dimension) {
	if (!in_domain(domain, val))
	    return NULL;
    } else {
	if (!in_domain(domain, val, dim, MaxRays))
	    return NULL;
    }
    Vector *res = Vector_Alloc(max.size());
    for (int i = 0; i < max.size(); ++i) {
	compute_evalue(max[i], val, &res->p[i]);
    }
    return res;
}

static Matrix *add_ge_constraint(EDomain *ED, evalue *constraint,
				 vector<evalue *>& new_floors)
{
    Polyhedron *D = ED->D;
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
    int fract = 0;
    for (evalue *e = constraint; value_zero_p(e->d); 
	 e = &e->x.p->arr[type_offset(e->x.p)]) {
	int i;
	if (e->x.p->type != fractional)
	    continue;
	for (i = 0; i < ED->floors.size(); ++i)
	    if (eequal(&e->x.p->arr[0], ED->floors[i]))
		break;
	if (i < ED->floors.size())
	    continue;
	++fract;
    }

    int rows = D->NbConstraints+2*fract+1;
    int cols = 2+D->Dimension+fract;
    Matrix *M = Matrix_Alloc(rows, cols);
    for (int i = 0; i < D->NbConstraints; ++i) {
	Vector_Copy(D->Constraint[i], M->p[i], 1+D->Dimension);
	value_assign(M->p[i][1+D->Dimension+fract], 
		     D->Constraint[i][1+D->Dimension]);
    }
    value_set_si(M->p[rows-1][0], 1);
    fract = 0;
    evalue *e;
    for (e = constraint; value_zero_p(e->d); e = &e->x.p->arr[type_offset(e->x.p)]) {
	if (e->x.p->type == fractional) {
	    int i, pos;

	    i = ED->find_floor(&e->x.p->arr[0]);
	    if (i >= 0)
		pos = D->Dimension-ED->floors.size()+i;
	    else
		pos = D->Dimension+fract;

	    add_fract(e, M->p[rows-1], cols, 1+pos);

	    if (pos < D->Dimension)
		continue;

	    /* constraints for the new floor */
	    int row = D->NbConstraints+2*fract;
	    value_set_si(M->p[row][0], 1);
	    evalue2constraint_r(NULL, &e->x.p->arr[0], M->p[row], cols);
	    value_oppose(M->p[row][1+D->Dimension+fract], M->p[row][0]);
	    value_set_si(M->p[row][0], 1);

	    Vector_Scale(M->p[row]+1, M->p[row+1]+1, mone.x.n, cols-1);
	    value_set_si(M->p[row+1][0], 1);
	    value_addto(M->p[row+1][cols-1], M->p[row+1][cols-1],
			M->p[row+1][1+D->Dimension+fract]);
	    value_decrement(M->p[row+1][cols-1], M->p[row+1][cols-1]);

	    evalue *arg = new evalue;
	    value_init(arg->d);
	    evalue_copy(arg, &e->x.p->arr[0]);
	    new_floors.push_back(arg);

	    ++fract;
	} else {
	    assert(e->x.p->type == polynomial);
	    assert(e->x.p->size == 2);
	    add_coeff(M->p[rows-1], cols, &e->x.p->arr[1], e->x.p->pos);
	}
    }
    add_coeff(M->p[rows-1], cols, e, cols-1);
    value_set_si(M->p[rows-1][0], 1);
    free_evalue_refs(&mone); 
    return M;
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

static Matrix *remove_equalities(Polyhedron **P, unsigned nparam, unsigned MaxRays);

Vector *Polyhedron_not_empty(Polyhedron *P, unsigned MaxRays)
{
    Polyhedron *Porig = P;
    Vector *sample;

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
		     unsigned MaxRays)
{
    Matrix *M, *M2;
    EDomain *EDlt = NULL, *EDeq = NULL, *EDgt = NULL;
    Polyhedron *D2;
    Value mone;
    value_init(mone);
    value_set_si(mone, -1);
    *Dlt = NULL;
    *Deq = NULL;
    *Dgt = NULL;
    vector<evalue *> new_floors;
    M = add_ge_constraint(D, sp.constraint, new_floors);
    if (sp.sign == split::lge || sp.sign == split::ge) {
	M2 = Matrix_Copy(M);
	value_decrement(M2->p[M2->NbRows-1][M2->NbColumns-1],
			M2->p[M2->NbRows-1][M2->NbColumns-1]);
	D2 = Constraints2Polyhedron(M2, MaxRays);
	EDgt = new EDomain(D2, D, new_floors);
	Polyhedron_Free(D2);
	Matrix_Free(M2);
    }
    if (sp.sign == split::lge || sp.sign == split::le) {
	M2 = Matrix_Copy(M);
	Vector_Scale(M2->p[M2->NbRows-1]+1, M2->p[M2->NbRows-1]+1,
		     mone, M2->NbColumns-1);
	value_decrement(M2->p[M2->NbRows-1][M2->NbColumns-1],
			M2->p[M2->NbRows-1][M2->NbColumns-1]);
	D2 = Constraints2Polyhedron(M2, MaxRays);
	EDlt = new EDomain(D2, D, new_floors);
	Polyhedron_Free(D2);
	Matrix_Free(M2);
    }

    assert(sp.sign == split::lge || sp.sign == split::ge || sp.sign == split::le);
    value_set_si(M->p[M->NbRows-1][0], 0);
    D2 = Constraints2Polyhedron(M, MaxRays);
    EDeq = new EDomain(D2, D, new_floors);
    Polyhedron_Free(D2);
    Matrix_Free(M);

    Vector *sample = D->sample;
    if (sample && new_floors.size() > 0) {
	assert(sample->Size == D->D->Dimension+1);
	sample = Vector_Alloc(D->D->Dimension+new_floors.size()+1);
	Vector_Copy(D->sample->p, sample->p, D->D->Dimension);
	value_set_si(sample->p[D->D->Dimension+new_floors.size()], 1);
	for (int i = 0; i < new_floors.size(); ++i)
	    compute_evalue(new_floors[i], sample->p, sample->p+D->D->Dimension+i);
    }

    for (int i = 0; i < new_floors.size(); ++i) {
	free_evalue_refs(new_floors[i]);
	delete new_floors[i];
    }

    if (EDeq) {
	if (sample && in_domain(EDeq->D, sample->p, sample->Size-1, MaxRays)) {
	    EDeq->sample = Vector_Alloc(sample->Size);
	    Vector_Copy(sample->p, EDeq->sample->p, sample->Size);
	} else if (!(EDeq->sample = Polyhedron_not_empty(EDeq->D, MaxRays))) {
	    delete EDeq;
	    EDeq = NULL;
	}
    }
    if (EDgt) {
	if (sample && in_domain(EDgt->D, sample->p, sample->Size-1, MaxRays)) {
	    EDgt->sample = Vector_Alloc(sample->Size);
	    Vector_Copy(sample->p, EDgt->sample->p, sample->Size);
	} else if (!(EDgt->sample = Polyhedron_not_empty(EDgt->D, MaxRays))) {
	    delete EDgt;
	    EDgt = NULL;
	}
    }
    if (EDlt) {
	if (sample && in_domain(EDlt->D, sample->p, sample->Size-1, MaxRays)) {
	    EDlt->sample = Vector_Alloc(sample->Size);
	    Vector_Copy(sample->p, EDlt->sample->p, sample->Size);
	} else if (!(EDlt->sample = Polyhedron_not_empty(EDlt->D, MaxRays))) {
	    delete EDlt;
	    EDlt = NULL;
	}
    }
    *Dlt = EDlt;
    *Deq = EDeq;
    *Dgt = EDgt;
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

static vector<max_term*> lexmin(indicator& ind, EDomain *D, unsigned nparam,
				unsigned MaxRays, vector<int> loc)
{
    vector<max_term*> maxima;
    int len = 1 + D->D->Dimension + 1;
    Value lcm, a, b;
    evalue mone;
    EDomain *Dorig = D;

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
    value_init(lcm);
    value_init(a);
    value_init(b);
    Vector *c = Vector_Alloc(len);
    Matrix *T = Matrix_Alloc(2, len-1);
    for (int i = 0; i < ind.term.size(); ++i) {
	bool restart = false;	/* true if we have modified ind from i up */
	bool stop = false;  	/* true if i can never be smallest */
	int peel = -1;	    	/* term to peel against */
	vector<split> splits;
	if (ind.term[i]->sign < 0)
	    continue;
	int dim = ind.term[i]->den.NumCols();
	int j;
	for (j = 0; j < ind.term.size(); ++j) {
	    if (i == j)
		continue;
	    int k;
	    for (k = 0; k < dim; ++k) {
		/* compute ind.term->[i]->vertex[k] - ind.term->[j]->vertex[k] */
		evalue *diff = new evalue;
		value_init(diff->d);
		evalue_copy(diff, ind.term[j]->vertex[k]);
		emul(&mone, diff);
		eadd(ind.term[i]->vertex[k], diff);
		reduce_evalue(diff);
		int fract = evalue2constraint(D, diff, c->p, len);
		Vector_Copy(c->p+1, T->p[0], len-1);
		value_assign(T->p[1][len-2], c->p[0]);

		int min, max;
		interval_minmax(D->D, T, &min, &max, MaxRays);
		if (max < 0) {
		    free_evalue_refs(diff); 
		    delete diff;
		    break;
		}
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
		if (min > 0) {
		    free_evalue_refs(diff); 
		    delete diff;
		    stop = true;
		    break;
		}
		if (max == 0 && min == 0) {
		    if (!EVALUE_IS_ZERO(*diff)) {
			ind.substitute(diff);
			ind.simplify();
			restart = true;
		    }
		    free_evalue_refs(diff); 
		    delete diff;
		    if (restart)
			break;
		    continue;
		}
		if (min < 0 && max == 0)
		    splits.push_back(split(diff, split::le));
		else if (max > 0 && min == 0)
		    splits.push_back(split(diff, split::ge));
		else
		    splits.push_back(split(diff, split::lge));
		break;
	    }
	    if (k == dim && ind.term[j]->sign < 0)
		peel = j;
	    if (stop || restart)
		break;
	}
	if (restart) {
	    /* The ith entry may have been removed, so we have to consider
	     * it again.
	     */
	    --i;
	    for (j = 0; j < splits.size(); ++j) {
		free_evalue_refs(splits[j].constraint);
		delete splits[j].constraint;
	    }
	    continue;
	}
	if (stop) {
	    for (j = 0; j < splits.size(); ++j) {
		free_evalue_refs(splits[j].constraint);
		delete splits[j].constraint;
	    }
	    continue;
	}
	if (peel != -1) {
	    // ind.peel(i, peel);
	    ind.combine(i, peel);
	    ind.simplify();
	    i = -1;		    /* start over */
	    for (j = 0; j < splits.size(); ++j) {
		free_evalue_refs(splits[j].constraint);
		delete splits[j].constraint;
	    }
	    continue;
	} 
	if (splits.size() != 0) {
	    for (j = 0; j < splits.size(); ++j)
		if (splits[j].sign == split::le)
		    break;
	    if (j == splits.size())
		j = 0;
		EDomain *Dlt, *Deq, *Dgt;
		split_on(splits[j], D, &Dlt, &Deq, &Dgt, MaxRays);
		assert(Dlt || Deq || Dgt);
		if (Deq) {
		    loc.push_back(0);
		    indicator indeq(ind);
		    indeq.substitute(splits[j].constraint);
		    Polyhedron *P = Polyhedron_Project_Initial(Deq->D, nparam);
		    P = DomainConstraintSimplify(P, MaxRays);
		    indeq.reduce_in_domain(P);
		    Polyhedron_Free(P);
		    indeq.simplify();
		    vector<max_term*> maxeq = lexmin(indeq, Deq, nparam,
						     MaxRays, loc);
		    maxima.insert(maxima.end(), maxeq.begin(), maxeq.end());
		    loc.pop_back();
		    delete Deq;
		}
		if (Dgt) {
		    loc.push_back(1);
		    indicator indgt(ind);
		    Polyhedron *P = Polyhedron_Project_Initial(Dgt->D, nparam);
		    P = DomainConstraintSimplify(P, MaxRays);
		    indgt.reduce_in_domain(P);
		    Polyhedron_Free(P);
		    indgt.simplify();
		    vector<max_term*> maxeq = lexmin(indgt, Dgt, nparam,
						     MaxRays, loc);
		    maxima.insert(maxima.end(), maxeq.begin(), maxeq.end());
		    loc.pop_back();
		    delete Dgt;
		}
		if (Dlt) {
		    loc.push_back(-1);
		    Polyhedron *P = Polyhedron_Project_Initial(Dlt->D, nparam);
		    P = DomainConstraintSimplify(P, MaxRays);
		    ind.reduce_in_domain(P);
		    Polyhedron_Free(P);
		    ind.simplify();
		    if (D != Dorig)
			delete D;
		    D = Dlt;
		    if (splits.size() > 1) {
			vector<max_term*> maxeq = lexmin(ind, Dlt, nparam,
							 MaxRays, loc);
			maxima.insert(maxima.end(), maxeq.begin(), maxeq.end());
			for (j = 0; j < splits.size(); ++j) {
			    free_evalue_refs(splits[j].constraint);
			    delete splits[j].constraint;
			}
			break;
		    }
		}
	    /* the vertex turned out not to be minimal */
	    for (j = 0; j < splits.size(); ++j) {
		free_evalue_refs(splits[j].constraint);
		delete splits[j].constraint;
	    }
	    if (!Dlt)
		break;
	}
	    max_term *maximum = new max_term;
	    maxima.push_back(maximum);
	    maximum->dim = nparam;
	    maximum->domain = Polyhedron_Copy(D->D);
	    for (int j = 0; j < dim; ++j) {
		evalue *E = new evalue;
		value_init(E->d);
		evalue_copy(E, ind.term[i]->vertex[j]);
		if (evalue_frac2floor_in_domain(E, D->D))
		    reduce_evalue(E);
		maximum->max.push_back(E);
	    }
	    break;
    }
    Matrix_Free(T);
    Vector_Free(c);
    value_clear(lcm);
    value_clear(a);
    value_clear(b);
    free_evalue_refs(&mone); 
    if (D != Dorig)
	delete D;
    return maxima;
}

static Matrix *compress_parameters(Polyhedron **P, unsigned nparam, unsigned MaxRays)
{
    Matrix *M, *T, *CP;

    /* compress_parms doesn't like equalities that only involve parameters */
    for (int i = 0; i < (*P)->NbEq; ++i)
	assert(First_Non_Zero((*P)->Constraint[i]+1, (*P)->Dimension-nparam) != -1);

    M = Matrix_Alloc((*P)->NbEq, (*P)->Dimension+2);
    Vector_Copy((*P)->Constraint[0], M->p[0], (*P)->NbEq * ((*P)->Dimension+2));
    CP = compress_parms(M, nparam);
    Matrix_Free(M);

    if (isIdentity(CP)) {
	Matrix_Free(CP);
	return NULL;
    }

    T = align_matrix(CP, (*P)->Dimension+1);
    *P = Polyhedron_Preimage(*P, T, MaxRays);
    Matrix_Free(T);

    return CP;
}

static Matrix *remove_equalities(Polyhedron **P, unsigned nparam, unsigned MaxRays)
{
    unsigned dim = (*P)->Dimension - nparam;
    Matrix *M, *H, *Q, *U, *C, *ratH, *invH, *Ul, *T1, *T2, *T;
    Value mone;
    int n;

    for (n = 0; n < (*P)->NbEq; ++n)
	if (First_Non_Zero((*P)->Constraint[n]+1, dim) == -1)
	    break;
    if (n == 0)
	return NULL;
    value_init(mone);
    value_set_si(mone, -1);
    M = Matrix_Alloc(n, dim);
    C = Matrix_Alloc(n+1, nparam+1);
    for (int i = 0; i < n; ++i) {
	Vector_Copy((*P)->Constraint[i]+1, M->p[i], dim);
	Vector_Scale((*P)->Constraint[i]+1+dim, C->p[i], mone, nparam+1);
    }
    value_set_si(C->p[n][nparam], 1);
    left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    value_clear(mone);

    /* we will need to treat scalings later */
    if (nparam > 0)
	for (int i = 0; i < n; ++i)
	    assert(value_one_p(H->p[i][i]));

    ratH = Matrix_Alloc(n+1, n+1);
    invH = Matrix_Alloc(n+1, n+1);
    for (int i = 0; i < n; ++i)
	Vector_Copy(H->p[i], ratH->p[i], n);
    value_set_si(ratH->p[n][n], 1);
    int ok = Matrix_Inverse(ratH, invH);
    Matrix_Free(H);
    Matrix_Free(ratH);
    assert(ok);
    T1 = Matrix_Alloc(n+1, nparam+1);
    Matrix_Product(invH, C, T1);
    if (nparam == 0 && value_notone_p(T1->p[n][nparam])) {
	for (int i = 0; i < n; ++i) {
	    if (!mpz_divisible_p(T1->p[i][nparam], T1->p[n][nparam])) {
		Matrix_Free(T1);
		Matrix_Free(C);
		Matrix_Free(invH);
		Matrix_Free(U);
		Matrix_Free(Q);
		return NULL;
	    }
	    value_division(T1->p[i][nparam], T1->p[i][nparam], T1->p[n][nparam]);
	}
	value_set_si(T1->p[n][nparam], 1);
    }
    Matrix_Free(C);
    Matrix_Free(invH);
    Ul = Matrix_Alloc(dim+1, n+1);
    for (int i = 0; i < dim; ++i)
	Vector_Copy(U->p[i], Ul->p[i], n);
    value_set_si(Ul->p[dim][n], 1);
    T2 = Matrix_Alloc(dim+1, nparam+1);
    Matrix_Product(Ul, T1, T2);
    Matrix_Free(Ul);
    Matrix_Free(T1);

    T = Matrix_Alloc(dim+1, (dim-n)+nparam+1);
    for (int i = 0; i < dim; ++i) {
	Vector_Copy(U->p[i]+n, T->p[i], dim-n);
	Vector_Copy(T2->p[i], T->p[i]+dim-n, nparam+1);
    }
    value_set_si(T->p[dim][(dim-n)+nparam], 1);
    assert(value_one_p(T2->p[dim][nparam]));
    Matrix_Free(U);
    Matrix_Free(T2);

    T2 = Matrix_Alloc((dim-n)+nparam+1, dim+nparam+1);
    for (int i = 0; i < dim-n; ++i)
	Vector_Copy(Q->p[n+i], T2->p[i], dim);
    for (int i = 0; i < nparam+1; ++i)
	value_set_si(T2->p[dim-n+i][dim+i], 1);
    *P = Polyhedron_Image(*P, T2, MaxRays);
    Matrix_Free(Q);
    Matrix_Free(T2);

    return T;
}

static vector<max_term*> lexmin(Polyhedron *P, Polyhedron *C, unsigned MaxRays)
{
    unsigned nparam = C->Dimension;
    Param_Polyhedron *PP = NULL;
    Polyhedron *CEq = NULL, *pVD;
    Matrix *CT = NULL;
    Matrix *T = NULL, *CP = NULL;
    Param_Domain *D, *next;
    Param_Vertices *V;
    Polyhedron *Porig = P;
    int i;
    vector<max_term*> all_max;
    Polyhedron *Q;

    if (emptyQ2(P))
	return all_max;

    POL_ENSURE_VERTICES(P);

    if (emptyQ2(P))
	return all_max;

    if (P->NbEq > 0) {
	if (nparam > 0)
	    CP = compress_parameters(&P, nparam, MaxRays);
	Q = P;
	T = remove_equalities(&P, nparam, MaxRays);
	if (Q != Porig)
	    Polyhedron_Free(Q);
    }

    Q = P;
    PP = Polyhedron2Param_SimplifiedDomain(&P,C,
					   (MaxRays & POL_NO_DUAL) ? 0 : MaxRays,
					   &CEq,&CT);
    if (P != Q && Q != Porig)
	Polyhedron_Free(Q);

    if (CT) {
	if (isIdentity(CT)) {
	    Matrix_Free(CT);
	    CT = NULL;
	} else
	    nparam = CT->NbRows - 1;
    }

    unsigned dim = P->Dimension - nparam;

    int nd;
    for (nd = 0, D=PP->D; D; ++nd, D=D->next);
    Polyhedron **fVD = new Polyhedron*[nd];

    indicator_constructor ic(P, dim, PP->nbV);

    for (i = 0, V = PP->V; V; V = V->next, i++) {
	ic.decompose_at_vertex(V, i, MaxRays);
    }
    if (T) {
	ic.substitute(T);
	ic.normalize();
    }

    for (nd = 0, D=PP->D; D; D=next) {
	next = D->next;

	Polyhedron *rVD = reduce_domain(D->Domain, CT, CEq,
					fVD, nd, MaxRays);
	if (!rVD)
	    continue;

	pVD = CT ? DomainImage(rVD,CT,MaxRays) : rVD;

	indicator ind;

	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    for (int j = 0; j < ic.terms[_i].size(); ++j) {
		indicator_term *term = new indicator_term(*ic.terms[_i][j]);
		term->reduce_in_domain(pVD);
		ind.term.push_back(term);
	    }
	END_FORALL_PVertex_in_ParamPolyhedron;

	ind.simplify();

	EDomain epVD(pVD);
	vector<int> loc;
	vector<max_term*> maxima = lexmin(ind, &epVD, nparam, MaxRays, loc);
	if (CP)
	    for (int j = 0; j < maxima.size(); ++j)
		maxima[j]->substitute(CP, MaxRays);
	all_max.insert(all_max.end(), maxima.begin(), maxima.end());

	++nd;
	if (rVD != pVD)
	    Domain_Free(pVD);
	Domain_Free(rVD);
    }
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

    while ((c = getopt_long(argc, argv, "TAm:M:r:V", options, &ind)) != -1) {
	switch (c) {
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
    C = Constraints2Polyhedron(MA, MAXRAYS);
    Matrix_Free(MA);
    fscanf(stdin, " %d", &bignum);
    assert(bignum == -1);
    MA = Matrix_Read();
    A = Constraints2Polyhedron(MA, MAXRAYS);
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

    iter_names = util_generate_names(A->Dimension - C->Dimension, "i");
    param_names = util_generate_names(C->Dimension, "p");
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    vector<max_term*> maxima = lexmin(A, C, MAXRAYS);
    for (int i = 0; i < maxima.size(); ++i)
	maxima[i]->print(cout, param_names);

    if (verify)
	verify_results(A, C, maxima, m, M, print_all, MAXRAYS);

    for (int i = 0; i < maxima.size(); ++i)
	delete maxima[i];

    util_free_names(A->Dimension - C->Dimension, iter_names);
    util_free_names(C->Dimension, param_names);
    Polyhedron_Free(A);
    Polyhedron_Free(C);

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
    CC2 = DomainIntersection(C, CC, MAXRAYS);
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
	CC2 = AddConstraints(MM->p[0], 2*CC->Dimension, CC, MAXRAYS);
	U = Universe_Polyhedron(0);
	CS = Polyhedron_Scan(CC2, U, MAXRAYS & POL_NO_DUAL ? 0 : MAXRAYS);
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

    S = Polyhedron_Scan(A, C, MAXRAYS & POL_NO_DUAL ? 0 : MAXRAYS);

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
