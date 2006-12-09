#define GMPRATIONAL
#include <setoper.h>
#include <cdd.h>
#include <barvinok/util.h>
#include "polysign.h"
#include "initcdd.h"

#define EMPTY_DOMAIN	-2

static int polyhedron_affine_minmax(DD_LPObjectiveType obj, Polyhedron *P,
				    Matrix *T, bool rational)
{
    DD_LPType  *lp;
    assert(P->Dimension == T->NbColumns-1);
    assert(T->NbRows == 2);
    DD_rowrange irev = P->NbConstraints;
    DD_rowrange rows = irev + P->NbEq + 1;
    DD_colrange cols = 1 + P->Dimension;
    lp = DD_CreateLPData(obj, DD_Rational, rows, cols);
    lp->Homogeneous = DD_FALSE;
    lp->objective = obj;

    for (DD_rowrange j = 0; j < P->NbConstraints; ++j) {
	for (DD_colrange k = 0; k < P->Dimension; ++k) {
	    DD_set_si(lp->A[j][1+k], 5);
	    DD_set_z(lp->A[j][1+k], P->Constraint[j][1+k]);
	}
	DD_set_z(lp->A[j][0], P->Constraint[j][1+P->Dimension]);
	if (j < P->NbEq) {
	    set_addelem(lp->equalityset, j+1);
	    for (DD_colrange k = 0; k < P->Dimension; ++k)
		DD_neg(lp->A[irev][1+k], lp->A[j][1+k]);
	    DD_neg(lp->A[irev][0], lp->A[j][0]);
	    ++irev;
	}
    }
    /* objective function */
    for (DD_colrange k = 0; k < P->Dimension; ++k)
	DD_set_z(lp->A[rows-1][1+k], T->p[0][k]);
    DD_set_z(lp->A[rows-1][0], T->p[0][P->Dimension]); 

    DD_ErrorType err = DD_NoError;
    DD_LPSolve(lp, DD_DualSimplex, &err);
    assert(err == DD_NoError);

    int sign;
    if (lp->LPS == DD_Optimal) {
	if (rational)
	    DD_rat_sign(sign, obj, lp->optvalue);
	else
	    /* The objective function has integer coefficients,
	     * so the optimal should be an integer (over the integer points)
	     */
	    DD_int_sign(sign, obj, lp->optvalue);
    } else if (lp->LPS == DD_DualInconsistent) {
	if (obj == DD_LPmin)
	    sign = -1;
	else
	    sign = 1;
    } else if (lp->LPS == DD_Inconsistent) {
	sign = EMPTY_DOMAIN;
    } else
	assert(0);

    DD_FreeLPData(lp);
    return sign;
}

enum order_sign cdd_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    if (emptyQ2(D))
	return order_undefined;

    INIT_CDD;
    bool rational = !POL_ISSET(options->MaxRays, POL_INTEGER);
    int min = polyhedron_affine_minmax(DD_LPmin, D, T, rational);
    if (min == EMPTY_DOMAIN)
	return order_undefined;
    if (min > 0)
	return order_gt;
    int max = polyhedron_affine_minmax(DD_LPmax, D, T, rational);
    assert(max != EMPTY_DOMAIN);
    if (max < 0)
	return order_lt;
    if (min == max)
	return order_eq;
    if (max == 0)
	return order_le;
    if (min == 0)
	return order_ge;
    return order_unknown;
}
