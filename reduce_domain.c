#include "reduce_domain.h"

Polyhedron *reduce_domain(Polyhedron *D, Matrix *CT, Polyhedron *CEq,
			  Polyhedron **fVD, int nd, unsigned MaxRays)
{
    Polyhedron *Dt, *rVD;
    int i;

    Dt = CT ? DomainPreimage(D, CT, MaxRays) : D;
    rVD = CEq ? DomainIntersection(Dt, CEq, MaxRays) : Domain_Copy(Dt);

    /* if rVD is empty or too small in geometric dimension */
    if(!rVD || emptyQ(rVD) ||
	    (CEq && rVD->Dimension-rVD->NbEq < Dt->Dimension-Dt->NbEq-CEq->NbEq)) {
	if(rVD)
	    Domain_Free(rVD);
	if (CT)
	    Domain_Free(Dt);
	return 0;		/* empty validity domain */
    }

    if (CT)
	Domain_Free(Dt);

    fVD[nd] = Domain_Copy(rVD);
    for (i = 0 ; i < nd; ++i) {
	Polyhedron *F;
	Polyhedron *I = DomainIntersection(fVD[nd], fVD[i], MaxRays);
	if (emptyQ(I)) {
	    Domain_Free(I);
	    continue;
	}
	F = DomainSimplify(I, fVD[nd], MaxRays);
	if (F->NbEq == 1) {
	    Polyhedron *T = rVD;
	    rVD = DomainDifference(rVD, F, MaxRays);
	    Domain_Free(T);
	}
	Domain_Free(F);
	Domain_Free(I);
    }

    rVD = DomainConstraintSimplify(rVD, MaxRays);
    if (emptyQ(rVD)) {
	Domain_Free(fVD[nd]);
	Domain_Free(rVD);
	return 0;
    }

    Value c;
    value_init(c);
    barvinok_count(rVD, &c, MaxRays);
    if (value_zero_p(c)) {
	Domain_Free(rVD);
	rVD = 0;
    }
    value_clear(c);

    return rVD;
}
