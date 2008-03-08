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

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

struct range_data {
    struct barvinok_options *options;
    unsigned		 nparam;
    const evalue	*poly;
    int 		*signs;
    int			 sign;
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

static void propagate_on_domain(Polyhedron *D, const evalue *poly,
				struct range_data *data)
{
    const evalue *save_poly = data->poly;

    assert(D->Dimension > data->nparam);

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
    return data.pl_all;
}
