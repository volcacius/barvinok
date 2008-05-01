#include <barvinok/options.h>
#include "bernoulli.h"
#include "euler.h"
#include "laurent.h"
#include "summate.h"
#include "section_array.h"

extern evalue *evalue_outer_floor(evalue *e);
extern int evalue_replace_floor(evalue *e, const evalue *floor, int var);
extern void evalue_drop_floor(evalue *e, const evalue *floor);

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static evalue *sum_over_polytope_with_equalities(Polyhedron *P, evalue *E,
				 unsigned nvar,
				 struct evalue_section_array *sections,
				 struct barvinok_options *options)
{
    unsigned dim = P->Dimension;
    unsigned new_dim, new_nparam;
    Matrix *T = NULL, *CP = NULL;
    evalue **subs;
    evalue *sum;
    int j;

    if (emptyQ(P))
	return evalue_zero();

    assert(P->NbEq > 0);

    remove_all_equalities(&P, NULL, &CP, &T, dim-nvar, options->MaxRays);

    if (emptyQ(P)) {
	Polyhedron_Free(P);
	return evalue_zero();
    }

    new_nparam = CP ? CP->NbColumns-1 : dim - nvar;
    new_dim = T ? T->NbColumns-1 : nvar + new_nparam;

    /* We can avoid these substitutions if E is a constant */
    subs = ALLOCN(evalue *, dim);
    for (j = 0; j < nvar; ++j) {
	if (T)
	    subs[j] = affine2evalue(T->p[j], T->p[nvar+new_nparam][new_dim],
				    new_dim);
	else
	    subs[j] = evalue_var(j);
    }
    for (j = 0; j < dim-nvar; ++j) {
	if (CP)
	    subs[nvar+j] = affine2evalue(CP->p[j], CP->p[dim-nvar][new_nparam],
					 new_nparam);
	else
	    subs[nvar+j] = evalue_var(j);
	evalue_shift_variables(subs[nvar+j], 0, new_dim-new_nparam);
    }

    E = evalue_dup(E);
    evalue_substitute(E, subs);
    reduce_evalue(E);

    for (j = 0; j < dim; ++j)
	evalue_free(subs[j]);
    free(subs);

    if (new_dim-new_nparam > 0) {
	sum = barvinok_sum_over_polytope(P, E, new_dim-new_nparam,
					 sections, options);
	evalue_free(E);
	Polyhedron_Free(P);
    } else {
	sum = ALLOC(evalue);
	value_init(sum->d);
	sum->x.p = new_enode(partition, 2, new_dim);
	EVALUE_SET_DOMAIN(sum->x.p->arr[0], P);
	value_clear(sum->x.p->arr[1].d);
	sum->x.p->arr[1] = *E;
	free(E);
    }

    if (CP) {
	evalue_backsubstitute(sum, CP, options->MaxRays);
	Matrix_Free(CP);
    }

    if (T)
	Matrix_Free(T);

    return sum;
}

/* Add two constraints corresponding to floor = floor(e/d),
 *
 *	 e - d t       >= 0
 *	-e + d t + d-1 >= 0
 *
 * e is assumed to be an affine expression.
 */
Polyhedron *add_floor_var(Polyhedron *P, unsigned nvar, const evalue *floor,
				     struct barvinok_options *options)
{
    int i;
    unsigned dim = P->Dimension+1;
    Matrix *M = Matrix_Alloc(P->NbConstraints+2, 2+dim);
    Polyhedron *CP;
    Value *d = &M->p[0][1+nvar];
    evalue_extract_affine(floor, M->p[0]+1, M->p[0]+1+dim, d);
    value_oppose(*d, *d);
    value_set_si(M->p[0][0], 1);
    value_set_si(M->p[1][0], 1);
    Vector_Oppose(M->p[0]+1, M->p[1]+1, M->NbColumns-1);
    value_subtract(M->p[1][1+dim], M->p[1][1+dim], *d);
    value_decrement(M->p[1][1+dim], M->p[1][1+dim]);

    for (i = 0; i < P->NbConstraints; ++i) {
	Vector_Copy(P->Constraint[i], M->p[i+2], 1+nvar);
	Vector_Copy(P->Constraint[i]+1+nvar, M->p[i+2]+1+nvar+1, dim-nvar-1+1);
    }

    CP = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    return CP;
}

static evalue *evalue_add(evalue *a, evalue *b)
{
    if (!a)
	return b;
    if (!b)
	return a;
    eadd(a, b);
    evalue_free(a);
    return b;
}

/* Compute sum of a step-polynomial over a polytope by grouping
 * terms containing the same floor-expressions and introducing
 * new variables for each such expression.
 * In particular, while there is any floor-expression left,
 * the step-polynomial is split into a polynomial containing
 * the expression, which is then converted to a new variable,
 * and a polynomial not containing the expression.
 */
static evalue *sum_step_polynomial(Polyhedron *P, evalue *E, unsigned nvar,
				     struct barvinok_options *options)
{
    evalue *floor;
    evalue *cur = E;
    evalue *sum = NULL;
    evalue *t;

    while ((floor = evalue_outer_floor(cur))) {
	Polyhedron *CP;
	evalue *converted = evalue_dup(cur);
	evalue *converted_floor = evalue_dup(floor);

	evalue_shift_variables(converted, nvar, 1);
	evalue_shift_variables(converted_floor, nvar, 1);
	evalue_replace_floor(converted, converted_floor, nvar);
	CP = add_floor_var(P, nvar, converted_floor, options);
	evalue_free(converted_floor);
	t = sum_step_polynomial(CP, converted, nvar+1, options);
	evalue_free(converted);
	Polyhedron_Free(CP);
	sum = evalue_add(t, sum);

	if (cur == E)
	    cur = evalue_dup(cur);
	evalue_drop_floor(cur, floor);
	evalue_free(floor);
    }

    if (EVALUE_IS_ZERO(*cur))
	t = NULL;
    else if (options->summation == BV_SUM_EULER)
	t = euler_summate(P, cur, nvar, options);
    else if (options->summation == BV_SUM_LAURENT)
	t = laurent_summate(P, cur, nvar, options);
    else
	assert(0);

    if (E != cur)
	evalue_free(cur);

    return evalue_add(t, sum);
}

evalue *barvinok_sum_over_polytope(Polyhedron *P, evalue *E, unsigned nvar,
				     struct evalue_section_array *sections,
				     struct barvinok_options *options)
{
    if (P->NbEq)
	return sum_over_polytope_with_equalities(P, E, nvar, sections, options);

    if (options->summation == BV_SUM_BERNOULLI)
	return bernoulli_summate(P, E, nvar, sections, options);
    else if (options->summation == BV_SUM_BARVINOK)
	return box_summate(P, E, nvar, options->MaxRays);

    evalue_frac2floor2(E, 0);

    return sum_step_polynomial(P, E, nvar, options);
}

evalue *barvinok_summate(evalue *e, int nvar, struct barvinok_options *options)
{
    int i;
    struct evalue_section_array sections;
    evalue *sum;

    assert(nvar >= 0);
    if (nvar == 0 || EVALUE_IS_ZERO(*e))
	return evalue_dup(e);

    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);

    evalue_section_array_init(&sections);
    sum = evalue_zero();

    for (i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *D;
	for (D = EVALUE_DOMAIN(e->x.p->arr[2*i]); D; D = D->next) {
	    Polyhedron *next = D->next;
	    evalue *tmp;
	    D->next = NULL;

	    tmp = barvinok_sum_over_polytope(D, &e->x.p->arr[2*i+1], nvar,
					     &sections, options);
	    assert(tmp);
	    eadd(tmp, sum);
	    evalue_free(tmp);

	    D->next = next;
	}
    }

    free(sections.s);

    reduce_evalue(sum);
    return sum;
}

evalue *evalue_sum(evalue *E, int nvar, unsigned MaxRays)
{
    evalue *sum;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    sum = barvinok_summate(E, nvar, options);
    barvinok_options_free(options);
    return sum;
}

evalue *esum(evalue *e, int nvar)
{
    evalue *sum;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    sum = barvinok_summate(e, nvar, options);
    barvinok_options_free(options);
    return sum;
}

/* Turn unweighted counting problem into "weighted" counting problem
 * with weight equal to 1 and call barvinok_summate on this weighted problem.
 */
evalue *barvinok_summate_unweighted(Polyhedron *P, Polyhedron *C,
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
    sum = barvinok_summate(&e, P->Dimension - C->Dimension, options);
    free_evalue_refs(&e);
    return sum;
}
