#include <assert.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "bernoulli.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))
#define REALLOCN(ptr,type,n) (type*)realloc(ptr, (n) * sizeof(type))

static struct bernoulli_coef bernoulli_coef;
static struct poly_list bernoulli;
static struct poly_list faulhaber;

struct bernoulli_coef *bernoulli_coef_compute(int n)
{
    int i, j;
    Value factor, tmp;

    if (n < bernoulli_coef.n)
	return &bernoulli_coef;

    if (n >= bernoulli_coef.size) {
	int size = 3*(n + 5)/2;
	Vector *b;

	b = Vector_Alloc(size);
	if (bernoulli_coef.n) {
	    Vector_Copy(bernoulli_coef.num->p, b->p, bernoulli_coef.n);
	    Vector_Free(bernoulli_coef.num);
	}
	bernoulli_coef.num = b;
	b = Vector_Alloc(size);
	if (bernoulli_coef.n) {
	    Vector_Copy(bernoulli_coef.den->p, b->p, bernoulli_coef.n);
	    Vector_Free(bernoulli_coef.den);
	}
	bernoulli_coef.den = b;
	b = Vector_Alloc(size);
	if (bernoulli_coef.n) {
	    Vector_Copy(bernoulli_coef.lcm->p, b->p, bernoulli_coef.n);
	    Vector_Free(bernoulli_coef.lcm);
	}
	bernoulli_coef.lcm = b;

	bernoulli_coef.size = size;
    }
    value_init(factor);
    value_init(tmp);
    for (i = bernoulli_coef.n; i <= n; ++i) {
	if (i == 0) {
	    value_set_si(bernoulli_coef.num->p[0], 1);
	    value_set_si(bernoulli_coef.den->p[0], 1);
	    value_set_si(bernoulli_coef.lcm->p[0], 1);
	    continue;
	}
	value_set_si(bernoulli_coef.num->p[i], 0);
	value_set_si(factor, -(i+1));
	for (j = i-1; j >= 0; --j) {
	    mpz_mul_ui(factor, factor, j+1);
	    mpz_divexact_ui(factor, factor, i+1-j);
	    value_division(tmp, bernoulli_coef.lcm->p[i-1],
			   bernoulli_coef.den->p[j]);
	    value_multiply(tmp, tmp, bernoulli_coef.num->p[j]);
	    value_multiply(tmp, tmp, factor);
	    value_addto(bernoulli_coef.num->p[i], bernoulli_coef.num->p[i], tmp);
	}
	mpz_mul_ui(bernoulli_coef.den->p[i], bernoulli_coef.lcm->p[i-1], i+1);
	value_gcd(tmp, bernoulli_coef.num->p[i], bernoulli_coef.den->p[i]);
	if (value_notone_p(tmp)) {
	    value_division(bernoulli_coef.num->p[i],
			    bernoulli_coef.num->p[i], tmp);
	    value_division(bernoulli_coef.den->p[i],
			    bernoulli_coef.den->p[i], tmp);
	}
	value_lcm(bernoulli_coef.lcm->p[i],
		  bernoulli_coef.lcm->p[i-1], bernoulli_coef.den->p[i]);
    }
    bernoulli_coef.n = n+1;
    value_clear(factor);
    value_clear(tmp);

    return &bernoulli_coef;
}

/*
 * Compute either Bernoulli B_n or Faulhaber F_n polynomials.
 *
 * B_n =         sum_{k=0}^n {  n  \choose k } b_k x^{n-k}
 * F_n = 1/(n+1) sum_{k=0}^n { n+1 \choose k } b_k x^{n+1-k}
 */
static struct poly_list *bernoulli_faulhaber_compute(int n, struct poly_list *pl,
						     int faulhaber)
{
    int i, j;
    Value factor;
    struct bernoulli_coef *bc;

    if (n < pl->n)
	return pl;

    if (n >= pl->size) {
	int size = 3*(n + 5)/2;
	Vector **poly;

	poly = ALLOCN(Vector *, size);
	for (i = 0; i < pl->n; ++i)
	    poly[i] = pl->poly[i];
	free(pl->poly);
	pl->poly = poly;

	pl->size = size;
    }

    bc = bernoulli_coef_compute(n);

    value_init(factor);
    for (i = pl->n; i <= n; ++i) {
	pl->poly[i] = Vector_Alloc(i+faulhaber+2);
	value_assign(pl->poly[i]->p[i+faulhaber], bc->lcm->p[i]);
	if (faulhaber)
	    mpz_mul_ui(pl->poly[i]->p[i+2], bc->lcm->p[i], i+1);
	else
	    value_assign(pl->poly[i]->p[i+1], bc->lcm->p[i]);
	value_set_si(factor, 1);
	for (j = 1; j <= i; ++j) {
	    mpz_mul_ui(factor, factor, i+faulhaber+1-j);
	    mpz_divexact_ui(factor, factor, j);
	    value_division(pl->poly[i]->p[i+faulhaber-j],
			   bc->lcm->p[i], bc->den->p[j]);
	    value_multiply(pl->poly[i]->p[i+faulhaber-j],
			   pl->poly[i]->p[i+faulhaber-j], bc->num->p[j]);
	    value_multiply(pl->poly[i]->p[i+faulhaber-j],
			   pl->poly[i]->p[i+faulhaber-j], factor);
	}
	Vector_Normalize(pl->poly[i]->p, i+faulhaber+2);
    }
    value_clear(factor);
    pl->n = n+1;

    return pl;
}

struct poly_list *bernoulli_compute(int n)
{
    return bernoulli_faulhaber_compute(n, &bernoulli, 0);
}

struct poly_list *faulhaber_compute(int n)
{
    return bernoulli_faulhaber_compute(n, &faulhaber, 1);
}

/* shift variables in polynomial one down */
static void shift(evalue *e)
{
    int i;
    if (value_notzero_p(e->d))
	return;
    assert(e->x.p->type == polynomial);
    assert(e->x.p->pos > 1);
    e->x.p->pos--;
    for (i = 0; i < e->x.p->size; ++i)
	shift(&e->x.p->arr[i]);
}

static evalue *shifted_copy(evalue *src)
{
    evalue *e = ALLOC(evalue);
    value_init(e->d);
    evalue_copy(e, src);
    shift(e);
    return e;
}

static evalue *power_sums(struct poly_list *faulhaber, evalue *poly,
			  Vector *arg, Value denom, int neg, int alt_neg)
{
    int i;
    evalue *base = affine2evalue(arg->p, denom, arg->Size-1);
    evalue *sum = evalue_zero();

    for (i = 1; i < poly->x.p->size; ++i) {
	evalue *term = evalue_polynomial(faulhaber->poly[i], base);
	evalue *factor = shifted_copy(&poly->x.p->arr[i]);
	emul(factor, term);
	if (alt_neg && (i % 2))
	    evalue_negate(term);
	eadd(term, sum);
	evalue_free(factor);
	evalue_free(term);
    }
    if (neg)
	evalue_negate(sum);
    evalue_free(base);

    return sum;
}

struct Bernoulli_data {
    unsigned MaxRays;
    struct evalue_section *s;
    int size;
    int ns;
    evalue *e;
};

static void Bernoulli_init(unsigned n, void *cb_data)
{
    struct Bernoulli_data *data = (struct Bernoulli_data *)cb_data;
    int cases = 5;

    if (cases * n <= data->size)
	return;

    data->size = cases * (n + 16);
    data->s = REALLOCN(data->s, struct evalue_section, data->size);
}

static void Bernoulli_cb(Matrix *M, Value *lower, Value *upper, void *cb_data)
{
    struct Bernoulli_data *data = (struct Bernoulli_data *)cb_data;
    Matrix *M2;
    Polyhedron *T;
    evalue *factor = NULL;
    evalue *linear = NULL;
    int constant = 0;
    Value tmp;
    unsigned dim = M->NbColumns-2;
    Vector *row;
    int cases = 5;

    assert(lower);
    assert(upper);
    assert(data->ns + cases <= data->size);

    M2 = Matrix_Copy(M);
    T = Constraints2Polyhedron(M2, data->MaxRays);
    Matrix_Free(M2);

    POL_ENSURE_VERTICES(T);
    if (emptyQ(T)) {
	Polyhedron_Free(T);
	return;
    }

    assert(lower != upper);

    row = Vector_Alloc(dim+1);
    value_init(tmp);
    if (value_notzero_p(data->e->d)) {
	factor = data->e;
	constant = 1;
    } else {
	assert(data->e->x.p->type == polynomial);
	if (data->e->x.p->pos > 1) {
	    factor = shifted_copy(data->e);
	    constant = 1;
	} else
	    factor = shifted_copy(&data->e->x.p->arr[0]);
    }
    if (!EVALUE_IS_ZERO(*factor)) {
	value_absolute(tmp, upper[1]);
	/* upper - lower */
	Vector_Combine(lower+2, upper+2, row->p, tmp, lower[1], dim+1);
	value_multiply(tmp, tmp, lower[1]);
	/* upper - lower + 1 */
	value_addto(row->p[dim], row->p[dim], tmp);

	linear = affine2evalue(row->p, tmp, dim);
	emul(factor, linear);
    } else
	linear = evalue_zero();

    if (constant) {
	data->s[data->ns].E = linear;
	data->s[data->ns].D = T;
	++data->ns;
    } else {
	evalue *poly_u = NULL, *poly_l = NULL;
	Polyhedron *D;
	struct poly_list *faulhaber;
	assert(data->e->x.p->type == polynomial);
	assert(data->e->x.p->pos == 1);
	faulhaber = faulhaber_compute(data->e->x.p->size-1);
	/* lower is the constraint
	 *			 b i - l' >= 0		i >= l'/b = l
	 * upper is the constraint
	 *			-a i + u' >= 0		i <= u'/a = u
	 */
	M2 = Matrix_Alloc(3, 2+T->Dimension);
	value_set_si(M2->p[0][0], 1);
	value_set_si(M2->p[1][0], 1);
	value_set_si(M2->p[2][0], 1);
	/* Case 1:
	 *		1 <= l		l' - b >= 0
	 */
	Vector_Oppose(lower+2, M2->p[0]+1, T->Dimension+1);
	value_subtract(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension],
		      lower[1]);
	D = AddConstraints(M2->p_Init, 1, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    evalue *extra;
	    if (!poly_u) {
		Vector_Copy(upper+2, row->p, dim+1);
		value_oppose(tmp, upper[1]);
		value_addto(row->p[dim], row->p[dim], tmp);
		poly_u = power_sums(faulhaber, data->e, row, tmp, 0, 0);
	    }
	    Vector_Oppose(lower+2, row->p, dim+1);
	    extra = power_sums(faulhaber, data->e, row, lower[1], 1, 0);
	    eadd(poly_u, extra);
	    eadd(linear, extra);

	    data->s[data->ns].E = extra;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 2:
	 *		1 <= -u		-u' - a >= 0
	 */
	Vector_Oppose(upper+2, M2->p[0]+1, T->Dimension+1);
	value_addto(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension],
		      upper[1]);
	D = AddConstraints(M2->p_Init, 1, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    evalue *extra;
	    if (!poly_l) {
		Vector_Copy(lower+2, row->p, dim+1);
		value_addto(row->p[dim], row->p[dim], lower[1]);
		poly_l = power_sums(faulhaber, data->e, row, lower[1], 0, 1);
	    }
	    Vector_Oppose(upper+2, row->p, dim+1);
	    value_oppose(tmp, upper[1]);
	    extra = power_sums(faulhaber, data->e, row, tmp, 1, 1);
	    eadd(poly_l, extra);
	    eadd(linear, extra);

	    data->s[data->ns].E = extra;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 3:
	 *		u >= 0		u' >= 0
	 *		-l >= 0		-l' >= 0
	 */
	Vector_Copy(upper+2, M2->p[0]+1, T->Dimension+1);
	Vector_Copy(lower+2, M2->p[1]+1, T->Dimension+1);
	D = AddConstraints(M2->p_Init, 2, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    if (!poly_l) {
		Vector_Copy(lower+2, row->p, dim+1);
		value_addto(row->p[dim], row->p[dim], lower[1]);
		poly_l = power_sums(faulhaber, data->e, row, lower[1], 0, 1);
	    }
	    if (!poly_u) {
		Vector_Copy(upper+2, row->p, dim+1);
		value_oppose(tmp, upper[1]);
		value_addto(row->p[dim], row->p[dim], tmp);
		poly_u = power_sums(faulhaber, data->e, row, tmp, 0, 0);
	    }
	
	    data->s[data->ns].E = ALLOC(evalue);
	    value_init(data->s[data->ns].E->d);
	    evalue_copy(data->s[data->ns].E, poly_u);
	    eadd(poly_l, data->s[data->ns].E);
	    eadd(linear, data->s[data->ns].E);
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 4:
	 *		l < 1		-l' + b - 1 >= 0
	 *		0 < l		l' - 1 >= 0
	 */
	Vector_Copy(lower+2, M2->p[0]+1, T->Dimension+1);
	value_addto(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension], lower[1]);
	value_decrement(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension]);
	Vector_Oppose(lower+2, M2->p[1]+1, T->Dimension+1);
	value_decrement(M2->p[1][1+T->Dimension], M2->p[1][1+T->Dimension]);
	D = AddConstraints(M2->p_Init, 2, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    if (!poly_u) {
		Vector_Copy(upper+2, row->p, dim+1);
		value_oppose(tmp, upper[1]);
		value_addto(row->p[dim], row->p[dim], tmp);
		poly_u = power_sums(faulhaber, data->e, row, tmp, 0, 0);
	    }

	    eadd(linear, poly_u);
	    data->s[data->ns].E = poly_u;
	    poly_u = NULL;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 5:
	 * 		-u < 1		u' + a - 1 >= 0
	 *		0 < -u		-u' - 1 >= 0
	 *		l <= 0		-l' >= 0
	 */
	Vector_Copy(upper+2, M2->p[0]+1, T->Dimension+1);
	value_subtract(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension],
			upper[1]);
	value_decrement(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension]);
	Vector_Oppose(upper+2, M2->p[1]+1, T->Dimension+1);
	value_decrement(M2->p[1][1+T->Dimension], M2->p[1][1+T->Dimension]);
	Vector_Copy(lower+2, M2->p[2]+1, T->Dimension+1);
	D = AddConstraints(M2->p_Init, 3, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    if (!poly_l) {
		Vector_Copy(lower+2, row->p, dim+1);
		value_addto(row->p[dim], row->p[dim], lower[1]);
		poly_l = power_sums(faulhaber, data->e, row, lower[1], 0, 1);
	    }

	    eadd(linear, poly_l);
	    data->s[data->ns].E = poly_l;
	    poly_l = NULL;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	Matrix_Free(M2);
	Polyhedron_Free(T);
	if (poly_l)
	    evalue_free(poly_l);
	if (poly_u)
	    evalue_free(poly_u);
	evalue_free(linear);
    }
    if (factor != data->e)
	evalue_free(factor);
    value_clear(tmp);
    Vector_Free(row);
}

/* Looks for variable with integer bounds, i.e., with coefficients 0, 1 or -1.
 * Returns 1 if such a variable is found and puts it in the first position,
 * possibly changing *P_p and *E_p.
 */
static int find_integer_bounds(Polyhedron **P_p, evalue **E_p, unsigned nvar)
{
    Polyhedron *P = *P_p;
    evalue *E = *E_p;
    unsigned dim = P->Dimension;
    int i, j;

    for (i = 0; i < nvar; ++i) {
	for (j = 0; j < P->NbConstraints; ++j) {
	    if (value_zero_p(P->Constraint[j][1+i]))
		continue;
	    if (value_one_p(P->Constraint[j][1+i]))
		continue;
	    if (value_mone_p(P->Constraint[j][1+i]))
		continue;
	    break;
	}
	if (j == P->NbConstraints)
	    break;
    }
    if (i == nvar)
	return 0;
    if (i == 0)
	return 1;
    P = Polyhedron_Copy(P);
    Polyhedron_ExchangeColumns(P, 1, 1+i);
    *P_p = P;

    if (value_zero_p(E->d)) {
	evalue **subs;
	subs = ALLOCN(evalue *, dim);
	for (j = 0; j < dim; ++j)
	    subs[j] = evalue_var(j);
	E = subs[0];
	subs[0] = subs[i];
	subs[i] = E;
	E = evalue_dup(*E_p);
	evalue_substitute(E, subs);
	for (j = 0; j < dim; ++j)
	    evalue_free(subs[j]);
	free(subs);
	*E_p = E;
    }

    return 1;
}

static evalue *sum_over_polytope(Polyhedron *P, evalue *E, unsigned nvar,
				 struct Bernoulli_data *data,
				 struct barvinok_options *options)
{
    unsigned dim = P->Dimension - 1;
    evalue *res;

    if (value_zero_p(P->Constraint[0][0]) &&
	    value_notzero_p(P->Constraint[0][1])) {
	res = ALLOC(evalue);
	value_init(res->d);
	value_set_si(res->d, 0);
	res->x.p = new_enode(partition, 2, dim);
	EVALUE_SET_DOMAIN(res->x.p->arr[0], Polyhedron_Project(P, dim));
	evalue_copy(&res->x.p->arr[1], E);
	reduce_evalue_in_domain(&res->x.p->arr[1], P);
	shift(&res->x.p->arr[1]);
    } else {
	data->ns = 0;
	data->e = E;

	for_each_lower_upper_bound(P, Bernoulli_init, Bernoulli_cb, data);

	res = evalue_from_section_array(data->s, data->ns);
    }

    if (nvar > 1) {
	evalue *tmp = Bernoulli_sum_evalue(res, nvar-1, options);
	evalue_free(res);
	res = tmp;
    }

    return res;
}

evalue *Bernoulli_sum_evalue(evalue *e, unsigned nvar,
			     struct barvinok_options *options)
{
    struct Bernoulli_data data;
    int i, j;
    evalue *sum = evalue_zero();

    if (EVALUE_IS_ZERO(*e))
	return sum;

    if (nvar == 0) {
	eadd(e, sum);
	return sum;
    }

    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);

    data.size = 16;
    data.s = ALLOCN(struct evalue_section, data.size);
    data.MaxRays = options->MaxRays;

    for (i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *D;
	for (D = EVALUE_DOMAIN(e->x.p->arr[2*i]); D; D = D->next) {
	    evalue *E = &e->x.p->arr[2*i+1];
	    Polyhedron *P = D;
	    Polyhedron *next = D->next;
	    evalue *tmp;
	    int integer_bounds;

	    P->next = NULL;

	    integer_bounds = find_integer_bounds(&P, &E, nvar);
	    if (options->approximation_method == BV_APPROX_NONE &&
		!integer_bounds) {
		evalue_free(sum);
		sum = NULL;
	    } else {
		evalue *tmp = sum_over_polytope(P, E, nvar, &data, options);
		eadd(tmp, sum);
		evalue_free(tmp);
	    }

	    if (P != D)
		Polyhedron_Free(P);
	    if (E != &e->x.p->arr[2*i+1])
		evalue_free(E);

	    D->next = next;;

	    if (!sum)
		break;
	}

	if (!sum)
	    break;
    }

    free(data.s);

    if (sum)
	reduce_evalue(sum);
    return sum;
}

evalue *Bernoulli_sum(Polyhedron *P, Polyhedron *C,
			struct barvinok_options *options)
{
    Polyhedron *CA, *D;
    evalue e;
    evalue *sum;

    if (emptyQ(P) || emptyQ(C))
	return evalue_zero();

    CA = align_context(C, P->Dimension, options->MaxRays);
    D = DomainIntersection(P, CA, options->MaxRays);
    Domain_Free(CA);

    if (emptyQ(D)) {
	Domain_Free(D);
	return evalue_zero();
    }

    value_init(e.d);
    e.x.p = new_enode(partition, 2, P->Dimension);
    EVALUE_SET_DOMAIN(e.x.p->arr[0], D);
    evalue_set_si(&e.x.p->arr[1], 1, 1);
    sum = Bernoulli_sum_evalue(&e, P->Dimension - C->Dimension, options);
    free_evalue_refs(&e);
    return sum;
}
