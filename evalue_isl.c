#include <isl_set_polylib.h>
#include <isl_constraint.h>
#include <barvinok/evalue.h>

static __isl_give isl_qpolynomial *extract_base(__isl_take isl_dim *dim,
	const evalue *e)
{
	int i;
	isl_vec *v;
	isl_div *div;
	isl_qpolynomial *base, *c;
	unsigned nparam;

	if (!dim)
		return NULL;

	if (e->x.p->type == polynomial)
		return isl_qpolynomial_var(dim, isl_dim_param, e->x.p->pos - 1);

	nparam = isl_dim_size(dim, isl_dim_param);
	v = isl_vec_alloc(dim->ctx, 2 + nparam);
	if (!v)
		goto error;

	isl_seq_clr(v->el + 2, nparam);
	evalue_extract_affine(&e->x.p->arr[0], v->el + 2, &v->el[1], &v->el[0]);

	div = isl_div_alloc(isl_dim_copy(dim));
	isl_div_set_constant(div, v->el[1]);
	isl_div_set_denominator(div, v->el[0]);

	for (i = 0; i < nparam; ++i)
		isl_div_set_coefficient(div, isl_dim_param, i, v->el[2 + i]);

	base = isl_qpolynomial_div(div);

	if (e->x.p->type == fractional) {
		base = isl_qpolynomial_neg(base);

		c = isl_qpolynomial_rat_cst(isl_dim_copy(dim), v->el[1], v->el[0]);
		base = isl_qpolynomial_add(base, c);

		for (i = 0; i < nparam; ++i) {
			isl_qpolynomial *t;
			c = isl_qpolynomial_rat_cst(isl_dim_copy(dim),
							v->el[2 + i], v->el[0]);
			t = isl_qpolynomial_var(isl_dim_copy(dim),
						isl_dim_param, i);
			t = isl_qpolynomial_mul(c, t);
			base = isl_qpolynomial_add(base, t);
		}
	}

	isl_vec_free(v);
	isl_dim_free(dim);

	return base;
error:
	isl_dim_free(dim);
	return NULL;
}

static int type_offset(enode *p)
{
	return p->type == fractional ? 1 : 
	       p->type == flooring ? 1 : 0;
}

static __isl_give isl_qpolynomial *evalue2qp(__isl_take isl_dim *dim,
	const evalue *e)
{
	int i;
	isl_qpolynomial *qp;
	isl_qpolynomial *base;
	int offset;

	if (EVALUE_IS_NAN(*e))
		return isl_qpolynomial_infty(dim);
	if (value_notzero_p(e->d))
		return isl_qpolynomial_rat_cst(dim, e->x.n, e->d);

	offset = type_offset(e->x.p);

	assert(e->x.p->type == polynomial ||
	       e->x.p->type == flooring ||
	       e->x.p->type == fractional);
	assert(e->x.p->size >= 1 + offset);

	base = extract_base(isl_dim_copy(dim), e);
	qp = evalue2qp(isl_dim_copy(dim), &e->x.p->arr[e->x.p->size - 1]);

	for (i = e->x.p->size - 2; i >= offset; --i) {
		qp = isl_qpolynomial_mul(qp, isl_qpolynomial_copy(base));
		qp = isl_qpolynomial_add(qp,
				    evalue2qp(isl_dim_copy(dim), &e->x.p->arr[i]));
	}

	isl_qpolynomial_free(base);
	isl_dim_free(dim);

	return qp;
}

static __isl_give isl_pw_qpolynomial *guarded_evalue2pwqp(__isl_take isl_set *set,
	const evalue *e);

static __isl_give isl_pw_qpolynomial *relation2pwqp(__isl_take isl_set *set,
	const evalue *e)
{
	int i;
	isl_vec *vec;
	isl_dim *dim;
	unsigned nparam;
	isl_pw_qpolynomial *pwqp;
	struct isl_constraint *c;
	struct isl_basic_set *bset;
	struct isl_set *guard;
	const evalue *fract;

	if (!set || !e)
		goto error;

	if (e->x.p->size == 1) {
		dim = isl_set_get_dim(set);
		isl_set_free(set);
		return isl_pw_qpolynomial_zero(dim);
	}

	isl_assert(set->ctx, e->x.p->size > 0, goto error);
	isl_assert(set->ctx, e->x.p->size <= 3, goto error);
	isl_assert(set->ctx, value_zero_p(e->x.p->arr[0].d), goto error);
	isl_assert(set->ctx, e->x.p->arr[0].x.p->type == fractional, goto error);
	fract = &e->x.p->arr[0];

	nparam = isl_set_dim(set, isl_dim_param);
	vec = isl_vec_alloc(set->ctx, 2 + nparam + 1);
	if (!vec)
		goto error;

	isl_seq_clr(vec->el + 2, nparam);
	evalue_extract_affine(&fract->x.p->arr[0],
				vec->el + 2, &vec->el[1], &vec->el[0]);

	dim = isl_set_get_dim(set);
	dim = isl_dim_add(dim, isl_dim_set, 1);

	bset = isl_basic_set_universe(dim);
	c = isl_equality_alloc(isl_dim_copy(dim));
	isl_int_neg(vec->el[0], vec->el[0]);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, vec->el[0]);
	isl_constraint_set_constant(c, vec->el[1]);
	for (i = 0; i < nparam; ++i)
		isl_constraint_set_coefficient(c, isl_dim_param, i, vec->el[2+i]);
	bset = isl_basic_set_add_constraint(bset, c);
	bset = isl_basic_set_project_out(bset, isl_dim_set, 0, 1);
	guard = isl_set_from_basic_set(bset);
	isl_vec_free(vec);

	pwqp = guarded_evalue2pwqp(isl_set_intersect(isl_set_copy(set),
						     isl_set_copy(guard)),
				   &e->x.p->arr[1]);

	if (e->x.p->size == 3) {
		isl_pw_qpolynomial *pwqpc;
		guard = isl_set_complement(guard);
		pwqpc = guarded_evalue2pwqp(isl_set_intersect(isl_set_copy(set),
							      isl_set_copy(guard)),
					    &e->x.p->arr[2]);
		pwqp = isl_pw_qpolynomial_add_disjoint(pwqp, pwqpc);
	}

	isl_set_free(set);
	isl_set_free(guard);

	return pwqp;
error:
	isl_set_free(set);
	return NULL;
}

static __isl_give isl_pw_qpolynomial *guarded_evalue2pwqp(__isl_take isl_set *set,
	const evalue *e)
{
	isl_qpolynomial *qp;

	if (value_zero_p(e->d) && e->x.p->type == relation)
		return relation2pwqp(set, e);

	qp = evalue2qp(isl_set_get_dim(set), e);

	return isl_pw_qpolynomial_alloc(set, qp);
}

__isl_give isl_pw_qpolynomial *isl_pw_qpolynomial_from_evalue(__isl_take isl_dim *dim, const evalue *e)
{
	int i;
	isl_pw_qpolynomial *pwqp;

	if (!dim || !e)
		goto error;
	if (EVALUE_IS_ZERO(*e))
		return isl_pw_qpolynomial_zero(dim);

	if (value_notzero_p(e->d)) {
		isl_set *set = isl_set_universe(isl_dim_copy(dim));
		isl_qpolynomial *qp = isl_qpolynomial_rat_cst(dim, e->x.n, e->d);
		return isl_pw_qpolynomial_alloc(set, qp);
	}

	assert(!EVALUE_IS_NAN(*e));

	assert(e->x.p->type == partition);

	pwqp = isl_pw_qpolynomial_zero(isl_dim_copy(dim));

	for (i = 0; i < e->x.p->size/2; ++i) {
		Polyhedron *D = EVALUE_DOMAIN(e->x.p->arr[2 * i]);
		isl_set *set = isl_set_new_from_polylib(D, isl_dim_copy(dim));
		isl_pw_qpolynomial *pwqp_i;

		pwqp_i = guarded_evalue2pwqp(set, &e->x.p->arr[2 * i + 1]);

		pwqp = isl_pw_qpolynomial_add_disjoint(pwqp, pwqp_i);
	}

	isl_dim_free(dim);

	return pwqp;
error:
	isl_dim_free(dim);
	return NULL;
}
