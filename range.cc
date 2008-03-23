#include <assert.h>
#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/bernstein.h>
#include <barvinok/options.h>
#include <barvinok/polylib.h>
#include <barvinok/util.h>
#include <barvinok/evalue.h>
#include "range.h"

using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

using std::cerr;
using std::endl;

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

struct range_data {
    struct barvinok_options *options;
    unsigned		 nparam;
    const evalue	*poly;
    int 		*signs;
    int			 sign;
    int			 test_monotonicity;
    int			 monotonicity;
    piecewise_lst	*pl;
    piecewise_lst	*pl_all;
    exvector	 	 params;
};

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 :
	  p->type == flooring ? 1 :
	  p->type == relation ? 1 : 0;
}

/* Remove all terms from e that are not positive (sign > 0)
 * or not negative (sign < 0)
 */
static void evalue_fixed_sign_terms(evalue *e, int* signs, int sign)
{
    int i, offset;
    int sign_odd = sign;

    assert(!EVALUE_IS_NAN(*e));
    if (value_notzero_p(e->d)) {
	if (sign * value_sign(e->x.n) < 0) {
	    value_set_si(e->x.n, 0);
	    value_set_si(e->d, 1);
	}
	return;
    }

    if (e->x.p->type == relation) {
	for (i = e->x.p->size-1; i >= 1; --i)
	    evalue_fixed_sign_terms(&e->x.p->arr[i], signs, sign);
	return;
    }

    if (e->x.p->type == polynomial)
	sign_odd *= signs[e->x.p->pos-1];

    offset = type_offset(e->x.p);
    for (i = offset; i < e->x.p->size; ++i)
	evalue_fixed_sign_terms(&e->x.p->arr[i], signs,
				     (i - offset) % 2 ? sign_odd : sign);
}

static void propagate_on_domain(Polyhedron *D, const evalue *poly,
				struct range_data *data);

static evalue *bound2evalue(Value *bound, unsigned dim)
{
    Value denom;
    evalue *b;

    value_init(denom);

    if (!bound) {
	b = evalue_nan();
    } else if (value_pos_p(bound[1])) {
	value_assign(denom, bound[1]);
	b = affine2evalue(bound+2, denom, dim);
	evalue_negate(b);
    } else {
	value_oppose(denom, bound[1]);
	b = affine2evalue(bound+2, denom, dim);
    }

    value_clear(denom);

    return b;
}

static int evalue_is_constant(const evalue *e)
{
    return EVALUE_IS_NAN(*e) || value_notzero_p(e->d);
}

static void range_cb(Matrix *M, Value *lower, Value *upper, void *cb_data)
{
    struct range_data *data = (struct range_data *)cb_data;
    evalue *l, *u;
    unsigned dim = M->NbColumns-2;
    evalue *pos, *neg;
    evalue **subs;
    Matrix *M2;
    Polyhedron *T;

    M2 = Matrix_Copy(M);
    T = Constraints2Polyhedron(M2, data->options->MaxRays);
    Matrix_Free(M2);

    POL_ENSURE_VERTICES(T);
    if (emptyQ(T)) {
	Polyhedron_Free(T);
	return;
    }

    l = bound2evalue(lower, dim);
    u = bound2evalue(upper, dim);

    subs = ALLOCN(evalue *, 1+dim);
    for (int i = 0; i < dim; ++i)
	subs[1+i] = evalue_var(i);

    if (evalue_is_constant(data->poly)) {
	pos = evalue_dup(data->poly);
    } else if (data->monotonicity) {
	pos = evalue_dup(data->poly);
	if (data->monotonicity * data->sign > 0)
	    subs[0] = u;
	else
	    subs[0] = l;
	evalue_substitute(pos, subs);
    } else {
	pos = evalue_dup(data->poly);
	neg = evalue_dup(data->poly);

	evalue_fixed_sign_terms(pos, data->signs, data->sign);
	evalue_fixed_sign_terms(neg, data->signs, -data->sign);

	subs[0] = u;
	evalue_substitute(pos, subs);

	subs[0] = l;
	evalue_substitute(neg, subs);

	eadd(neg, pos);
	evalue_free(neg);
    }

    for (int i = 0; i < dim; ++i)
	evalue_free(subs[1+i]);
    free(subs);
    evalue_free(l);
    evalue_free(u);

    if (dim == data->nparam) {
	data->pl->add_guarded_lst(T, lst(evalue2ex(pos, data->params)));
    } else {
	propagate_on_domain(T, pos, data);
	Polyhedron_Free(T);
    }

    evalue_free(pos);
}

static int has_sign(Polyhedron *D, evalue *poly, int sign,
		    int *signs, struct barvinok_options *options)
{
    struct range_data data_m;
    int i;
    lst l;

    data_m.options = options;
    data_m.nparam = 0;
    data_m.test_monotonicity = 0;
    data_m.signs = signs;

    data_m.pl_all = NULL;
    data_m.sign = -sign;
    propagate_on_domain(D, poly, &data_m);

    assert(data_m.pl_all->list.size() == 1);
    assert(universeQ(data_m.pl_all->list[0].first));
    l = data_m.pl_all->list[0].second;

    for (i = 0; i < l.nops(); ++i) {
	if (is_exactly_a<fail>(l.op(i)))
	    break;
	if (sign * l.op(i) < 0)
	    break;
    }
    delete data_m.pl_all;

    return i == l.nops();
}

/* Returns 1 if poly is monotonically increasing,
 *        -1 if poly is monotonically decreasing,
 *	   0 if no conclusion.
 */
static int monotonicity(Polyhedron *D, const evalue *poly,
			struct range_data *data)
{
    evalue **subs;
    evalue *diff;
    Value one;
    int result = 0;

    value_init(one);
    value_set_si(one, 1);

    subs = ALLOCN(evalue *, D->Dimension);
    for (int i = 0; i < D->Dimension; ++i)
	subs[i] = evalue_var(i);

    evalue_add_constant(subs[0], one);

    diff = evalue_dup(poly);
    evalue_substitute(diff, subs);

    for (int i = 0; i < D->Dimension; ++i)
	evalue_free(subs[i]);
    free(subs);

    evalue_negate(diff);
    eadd(poly, diff);
    reduce_evalue(diff);
    evalue_negate(diff);

    value_clear(one);

    if (has_sign(D, diff, 1, data->signs, data->options))
	result = 1;
    else if (has_sign(D, diff, -1, data->signs, data->options))
	result = -1;

    evalue_free(diff);
    return result;
}

static void propagate_on_domain(Polyhedron *D, const evalue *poly,
				struct range_data *data)
{
    const evalue *save_poly = data->poly;
    int save_monotonicity = data->monotonicity;

    assert(D->Dimension > data->nparam);

    if (data->test_monotonicity)
	data->monotonicity = monotonicity(D, poly, data);
    else
	data->monotonicity = 0;

    if (D->Dimension == data->nparam+1)
	data->pl = new piecewise_lst(data->params);

    data->poly = poly;
    for_each_lower_upper_bound(D, NULL, range_cb, data);

    if (D->Dimension == data->nparam+1) {
	if (!data->pl_all)
	    data->pl_all = data->pl;
	else {
	    data->pl_all->combine(*data->pl);
	    delete data->pl;
	}
	data->pl = NULL;
    }

    data->poly = save_poly;
    data->monotonicity = save_monotonicity;
}

/*
 * Determines the sign of each variable in D.
 * Copied from evalue.c
 */
static void domain_signs(Polyhedron *D, int *signs)
{
    int j, k;

    POL_ENSURE_VERTICES(D);
    for (j = 0; j < D->Dimension; ++j) {
	signs[j] = 0;
	for (k = 0; k < D->NbRays; ++k) {
	    signs[j] = value_sign(D->Ray[k][1+j]);
	    if (signs[j])
		break;
	}
    }
}

static evalue *ex2evalue(const ex& poly, const exvector& params, int pos)
{
    if (pos >= params.size()) {
	evalue *c = ALLOC(evalue);
	value_init(c->d);
	value_init(c->x.n);
	assert(is_a<numeric>(poly));
	numeric2value(ex_to<numeric>(poly).numer(), c->x.n);
	numeric2value(ex_to<numeric>(poly).denom(), c->d);
	return c;
    }

    evalue *v = evalue_var(pos);
    evalue *sum = ex2evalue(poly.coeff(params[pos], poly.degree(params[pos])),
			    params, pos+1);
    for (int i = poly.degree(params[pos])-1; i >= 0; --i) {
	evalue *t = ex2evalue(poly.coeff(params[pos], i), params, pos+1);
	emul(v, sum);
	eadd(t, sum);
	evalue_free(t);
    }
    evalue_free(v);
    return sum;
}

/* Returns true is poly is no better than any of those from begin to end */
static int is_no_better(const ex& poly, const lst::const_iterator& begin,
			const lst::const_iterator& end, const exvector& params,
			int sign, Polyhedron *D,
			int *signs, barvinok_options *options)
{
    lst::const_iterator k;
    int no_better = 0;

    for (k = begin; k != end; ++k) {
	ex diff = *k - poly;
	diff = diff.expand();
	evalue *e = ex2evalue(diff, params, 0);
	no_better = has_sign(D, e, sign, signs, options);
	evalue_free(e);
	if (no_better)
	    break;
    }
    return no_better;
}

static void remove_redundants(piecewise_lst *pl, const exvector& params,
				barvinok_options *options)
{
    if (pl->sign == 0)
	return;

    int *signs = (int *)alloca(sizeof(int) * params.size());

    for (int i = 0; i < pl->list.size(); ++i) {
	Polyhedron *D = pl->list[i].first;
	lst todo = pl->list[i].second;
	lst newlist;
	lst::const_iterator j, k;

	domain_signs(D, signs);

	for (j = todo.begin(); j != todo.end(); ++j) {
	    k = j; ++k;
	    if (is_no_better(*j, k, todo.end(), params,
				pl->sign, D, signs, options))
		continue;

	    if (is_no_better(*j, newlist.begin(), newlist.end(),
				params, pl->sign, D, signs, options))
		continue;

	    newlist.append(*j);
	}
	pl->list[i].second = newlist;
    }
}

piecewise_lst *evalue_range_propagation(piecewise_lst *pl_all, const evalue *e,
				      const exvector& params,
				      barvinok_options *options)
{
    evalue *e2;
    struct range_data data;

    if (EVALUE_IS_ZERO(*e))
	return pl_all;

    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);
    assert(e->x.p->size >= 2);
    unsigned nparam = params.size();
    unsigned dim = EVALUE_DOMAIN(e->x.p->arr[0])->Dimension;
    unsigned nvars = dim - nparam;

    if (nvars == 0) {
	assert(0);
	return pl_all;
    }

    data.nparam = nparam;
    data.params = params;
    data.options = options;
    data.pl_all = pl_all;
    if (options->bernstein_optimize == BV_BERNSTEIN_MIN)
	data.sign = -1;
    else
	data.sign = 1;
    data.test_monotonicity = 1;

    e2 = evalue_dup(e);
    evalue_split_domains_into_orthants(e2, options->MaxRays);

    data.signs = (int *)alloca(sizeof(int) * dim);

    for (int i = 0; i < e2->x.p->size/2; ++i) {
	Polyhedron *D = EVALUE_DOMAIN(e2->x.p->arr[2*i]);
	for (Polyhedron *P = D; P; P = P->next) {
	    domain_signs(P, data.signs);
	    propagate_on_domain(P, &e2->x.p->arr[2*i+1], &data);
	}
    }
    evalue_free(e2);
    if (data.pl_all) {
	data.pl_all->sign = options->bernstein_optimize;
	remove_redundants(data.pl_all, params, options);
    }
    return data.pl_all;
}
