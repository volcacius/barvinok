#include <barvinok/options.h>
#include "reduce_domain.h"

Polyhedron *reduce_domain(Polyhedron *D, Matrix *CT, Polyhedron *CEq,
			  Polyhedron **fVD, int nd,
			  struct barvinok_options *options)
{
    Polyhedron *Dt, *rVD;
    Polyhedron *C;
    Value c;
    int i;

    C = D->next ? DomainConvex(D, options->MaxRays) : D;
    Dt = CT ? DomainPreimage(C, CT, options->MaxRays) : C;
    rVD = CEq ? DomainIntersection(Dt, CEq, options->MaxRays) : Domain_Copy(Dt);

    /* if rVD is empty or too small in geometric dimension */
    if(!rVD || emptyQ(rVD) ||
	    (CEq && rVD->Dimension-rVD->NbEq < Dt->Dimension-Dt->NbEq-CEq->NbEq)) {
	if(rVD)
	    Domain_Free(rVD);
	if (D->next)
	    Polyhedron_Free(C);
	if (CT)
	    Domain_Free(Dt);
	return 0;		/* empty validity domain */
    }

    if (D->next)
	Polyhedron_Free(C);
    if (CT)
	Domain_Free(Dt);

    fVD[nd] = Domain_Copy(rVD);
    for (i = 0 ; i < nd; ++i) {
	Polyhedron *F;
	Polyhedron *I = DomainIntersection(fVD[nd], fVD[i], options->MaxRays);
	if (emptyQ(I)) {
	    Domain_Free(I);
	    continue;
	}
	F = DomainSimplify(I, fVD[nd], options->MaxRays);
	if (F->NbEq == 1) {
	    Polyhedron *T = rVD;
	    rVD = DomainDifference(rVD, F, options->MaxRays);
	    Domain_Free(T);
	}
	Domain_Free(F);
	Domain_Free(I);
    }

    rVD = DomainConstraintSimplify(rVD, options->MaxRays);
    if (emptyQ(rVD)) {
	Domain_Free(fVD[nd]);
	Domain_Free(rVD);
	return 0;
    }

    value_init(c);
    barvinok_count_with_options(rVD, &c, options);
    if (value_zero_p(c)) {
	Domain_Free(rVD);
	rVD = 0;
    }
    value_clear(c);

    return rVD;
}
