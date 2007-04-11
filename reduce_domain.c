#include <barvinok/options.h>
#include "reduce_domain.h"

Vector *inner_point(Polyhedron *P)
{
    Vector *average;
    int nv = 0;
    Value fc, fv;
    unsigned dim = P->Dimension;
    int i;

    average = Vector_Alloc(1+dim+1);

    POL_ENSURE_VERTICES(P);
    for (i = 0; i < P->NbRays; ++i) {
	if (value_zero_p(P->Ray[i][1+dim]))
	    continue;
	++nv;
	if (nv == 1) {
	    Vector_Copy(P->Ray[i]+1, average->p+1, dim+1);
	    continue;
	}
	if (nv == 2) {
	    value_init(fc);
	    value_init(fv);
	}
	value_assign(fc, average->p[1+dim]);
	value_lcm(fc, P->Ray[i][1+dim], &average->p[1+dim]);
	value_division(fc, average->p[1+dim], fc);
	value_division(fv, average->p[1+dim], P->Ray[i][1+dim]);
	Vector_Combine(average->p+1, P->Ray[i]+1, average->p+1, fc, fv, dim);
    }
    if (nv > 1) {
	value_set_si(fc, nv);
	value_multiply(average->p[1+dim], average->p[1+dim], fc);
	Vector_Normalize(average->p+1, dim+1);
	value_clear(fc);
	value_clear(fv);
    }
    for (i = 0; i < P->NbRays; ++i) {
	if (value_notzero_p(P->Ray[i][1+dim]))
	    continue;
	Vector_Add(average->p+1, P->Ray[i]+1, average->p+1, dim);
    }
    return average;
}

int is_internal(Vector *point, Value *constraint)
{
    int p;
    unsigned dim = point->Size-2;

    Inner_Product(constraint+1, point->p+1, dim+1, &point->p[0]);
    if (value_notzero_p(point->p[0]))
	return value_pos_p(point->p[0]);

    p = First_Non_Zero(constraint+1, dim);
    return value_pos_p(constraint[1+p]);
}

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

    if (CT)
	Domain_Free(Dt);

    fVD[nd] = Domain_Copy(C);
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
	    Polyhedron *FE = CT ? DomainPreimage(F, CT, options->MaxRays) : F;
	    rVD = DomainDifference(rVD, FE, options->MaxRays);
	    if (CT)
		Domain_Free(FE);
	    Domain_Free(T);
	}
	Domain_Free(F);
	Domain_Free(I);
    }

    if (D->next)
	Polyhedron_Free(C);

    rVD = DomainConstraintSimplify(rVD, options->MaxRays);
    if (emptyQ(rVD)) {
	Domain_Free(fVD[nd]);
	Domain_Free(rVD);
	return 0;
    }

    value_init(c);
    barvinok_count_with_options(rVD, &c, options);
    if (value_zero_p(c)) {
	Domain_Free(fVD[nd]);
	Domain_Free(rVD);
	rVD = 0;
    }
    value_clear(c);

    return rVD;
}
