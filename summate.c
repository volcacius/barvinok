#include <barvinok/options.h>
#include "bernoulli.h"
#include "euler.h"
#include "laurent.h"
#include "summate.h"
#include "section_array.h"

static evalue *sum_over_polytope(Polyhedron *P, evalue *E, unsigned nvar,
				 struct evalue_section_array *sections,
				 struct barvinok_options *options)
{
    if (options->summation == BV_SUM_EULER)
	return euler_summate(P, E, nvar, options);
    else if (options->summation == BV_SUM_LAURENT)
	return laurent_summate(P, E, nvar, options);
    else if (options->summation == BV_SUM_BERNOULLI)
	return bernoulli_summate(P, E, nvar, sections, options);
    else
	return box_summate(P, E, nvar, options->MaxRays);
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

	    tmp = sum_over_polytope(D, &e->x.p->arr[2*i+1], nvar,
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
