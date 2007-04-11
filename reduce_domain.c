#include <barvinok/options.h>
#include <barvinok/util.h>
#include "reduce_domain.h"

Polyhedron *true_context(Polyhedron *P, Matrix *CT,
			 Polyhedron *C, unsigned MaxRays)
{
    unsigned nparam = CT ? CT->NbRows - 1 : C->Dimension;
    Polyhedron *tmp = Polyhedron_Project(P, nparam);
    Polyhedron *D = CT ? DomainPreimage(tmp, CT, MaxRays) : tmp;
    C = DomainIntersection(D, C, MaxRays);
    if (CT)
	Polyhedron_Free(D);
    Polyhedron_Free(tmp);
    return C;
}

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

Polyhedron *reduce_domain(Polyhedron *D, Matrix *CT, Polyhedron *CEq, int nd,
			  Vector *inner, struct barvinok_options *options)
{
    Polyhedron *Dt, *rVD;
    Polyhedron *C;
    Value c;
    int i;
    Matrix *constraints;
    int changed;

    C = D->next ? DomainConvex(D, options->MaxRays) : D;
    Dt = CT ? DomainPreimage(C, CT, options->MaxRays) : C;
    rVD = CEq ? DomainIntersection(Dt, CEq, options->MaxRays) : Domain_Copy(Dt);

    /* If there is only one chamber, then we don't need to take care
     * of possible overlaps.
     * Plus, this decomposition may be the result of a recursive call
     * and then some of the assumptions used in determining whether
     * the domain is too small in geometric dimension no longer apply.
     */
    if (nd == 1)
	goto done;

    /* if rVD is empty or too small in geometric dimension */
    if(!rVD || emptyQ(rVD) ||
	    (CEq && rVD->Dimension-rVD->NbEq < Dt->Dimension-Dt->NbEq-CEq->NbEq)) {
	if (rVD)
	    Domain_Free(rVD);
	rVD = NULL;		/* empty validity domain */
done:
	if (D->next)
	    Polyhedron_Free(C);
	if (CT)
	    Domain_Free(Dt);
	return rVD;
    }

    if (CT)
	Domain_Free(Dt);

    assert(rVD->Dimension == inner->Size-2);
    constraints = Polyhedron2Constraints(rVD);
    changed = 0;
    for (i = 0; i < constraints->NbRows; ++i) {
	if (!is_internal(inner, constraints->p[i])) {
	    value_decrement(constraints->p[i][1+rVD->Dimension],
			    constraints->p[i][1+rVD->Dimension]);
	    changed = 1;
	}
    }
    if (changed) {
	Polyhedron_Free(rVD);
	rVD = Constraints2Polyhedron(constraints, options->MaxRays);
    }
    Matrix_Free(constraints);

    if (D->next)
	Polyhedron_Free(C);

    rVD = DomainConstraintSimplify(rVD, options->MaxRays);
    POL_ENSURE_FACETS(rVD);
    if (emptyQ(rVD)) {
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
