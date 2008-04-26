#include <barvinok/options.h>
#include "bernoulli.h"
#include "euler.h"
#include "laurent.h"
#include "summate.h"

evalue *barvinok_summate(evalue *e, int nvar, struct barvinok_options *options)
{
    if (options->summation == BV_SUM_EULER)
	return euler_summate(e, nvar, options);
    else if (options->summation == BV_SUM_LAURENT)
	return laurent_summate(e, nvar, options);
    else if (options->summation == BV_SUM_BERNOULLI)
	return Bernoulli_sum_evalue(e, nvar, options);
    else
	return evalue_sum(e, nvar, options->MaxRays);
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
