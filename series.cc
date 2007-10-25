#include <assert.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "genfun_constructor.h"
#include "remove_equalities.h"

/* Check whether all rays point in the positive directions
 * for the parameters
 */
static bool Polyhedron_has_positive_rays(Polyhedron *P, unsigned nparam)
{
    int r;
    for (r = 0; r < P->NbRays; ++r)
	if (value_zero_p(P->Ray[r][P->Dimension+1])) {
	    int i;
	    for (i = P->Dimension - nparam; i < P->Dimension; ++i)
		if (value_neg_p(P->Ray[r][i+1]))
		    return false;
	}
    return true;
}

/*
 * remove equalities that require a "compression" of the parameters
 */
static Polyhedron *remove_more_equalities(Polyhedron *P, unsigned nparam,
					  Matrix **CP, unsigned MaxRays)
{
    Polyhedron *Q = P;
    remove_all_equalities(&P, NULL, CP, NULL, nparam, MaxRays);
    if (P != Q)
	Polyhedron_Free(Q);
    return P;
}

/* frees P */
static gen_fun *series(Polyhedron *P, unsigned nparam, barvinok_options *options)
{
    Matrix *CP = NULL;
    gen_fun *gf;

    if (emptyQ2(P)) {
	Polyhedron_Free(P);
	return new gen_fun(Empty_Polyhedron(nparam));
    }

    assert(!Polyhedron_is_unbounded(P, nparam, options->MaxRays));
    assert(P->NbBid == 0);
    assert(Polyhedron_has_revlex_positive_rays(P, nparam));
    if (P->NbEq != 0)
	P = remove_more_equalities(P, nparam, &CP, options->MaxRays);
    assert(emptyQ2(P) || P->NbEq == 0);
    if (CP)
	nparam = CP->NbColumns-1;

    if (nparam == 0) {
	Value c;
	value_init(c);
	barvinok_count_with_options(P, &c, options);
	gf = new gen_fun(c);
	value_clear(c);
    } else {
	POL_ENSURE_VERTICES(P);
	if (P->NbEq)
	    gf = series(Polyhedron_Copy(P), nparam, options);
	else {
	    gf_base *red;
	    red = gf_base::create(Polyhedron_Project(P, nparam),
				  P->Dimension, nparam, options);
	    red->start_gf(P, options);
	    gf = red->gf;
	    delete red;
	}
    }
    if (CP) {
	gf->substitute(CP);
	Matrix_Free(CP);
    }
    Polyhedron_Free(P);
    return gf;
}

gen_fun * barvinok_series_with_options(Polyhedron *P, Polyhedron* C,
				       barvinok_options *options)
{
    Polyhedron *CA;
    unsigned nparam = C->Dimension;
    gen_fun *gf;

    CA = align_context(C, P->Dimension, options->MaxRays);
    P = DomainIntersection(P, CA, options->MaxRays);
    Polyhedron_Free(CA);

    gf = series(P, nparam, options);

    return gf;
}

gen_fun * barvinok_series(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    gen_fun *gf;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    gf = barvinok_series_with_options(P, C, options);
    barvinok_options_free(options);
    return gf;
}

static Polyhedron *skew_into_positive_orthant(Polyhedron *D, unsigned nparam, 
					      unsigned MaxRays)
{
    Matrix *M = NULL;
    Value tmp;
    value_init(tmp);
    for (Polyhedron *P = D; P; P = P->next) {
	POL_ENSURE_VERTICES(P);
	assert(!Polyhedron_is_unbounded(P, nparam, MaxRays));
	assert(P->NbBid == 0);
	assert(Polyhedron_has_positive_rays(P, nparam));

	for (int r = 0; r < P->NbRays; ++r) {
	    if (value_notzero_p(P->Ray[r][P->Dimension+1]))
		continue;
	    for (int i = 0; i < nparam; ++i) {
		int j;
		if (value_posz_p(P->Ray[r][i+1]))
		    continue;
		if (!M) {
		    M = Matrix_Alloc(D->Dimension+1, D->Dimension+1);
		    for (int i = 0; i < D->Dimension+1; ++i)
			value_set_si(M->p[i][i], 1);
		} else {
		    Inner_Product(P->Ray[r]+1, M->p[i], D->Dimension+1, &tmp);
		    if (value_posz_p(tmp))
			continue;
		}
		for (j = P->Dimension - nparam; j < P->Dimension; ++j)
		    if (value_pos_p(P->Ray[r][j+1]))
			break;
		assert(j < P->Dimension);
		value_pdivision(tmp, P->Ray[r][j+1], P->Ray[r][i+1]);
		value_subtract(M->p[i][j], M->p[i][j], tmp);
	    }
	}
    }
    value_clear(tmp);
    if (M) {
	D = DomainImage(D, M, MaxRays);
	Matrix_Free(M);
    }
    return D;
}

gen_fun* barvinok_enumerate_union_series_with_options(Polyhedron *D, Polyhedron* C, 
						      barvinok_options *options)
{
    Polyhedron *conv, *D2;
    Polyhedron *CA;
    gen_fun *gf = NULL, *gf2;
    unsigned nparam = C->Dimension;
    ZZ one, mone;
    one = 1;
    mone = -1;

    CA = align_context(C, D->Dimension, options->MaxRays);
    D = DomainIntersection(D, CA, options->MaxRays);
    Polyhedron_Free(CA);

    D2 = skew_into_positive_orthant(D, nparam, options->MaxRays);
    for (Polyhedron *P = D2; P; P = P->next) {
	assert(P->Dimension == D2->Dimension);
	gen_fun *P_gf;

	P_gf = series(Polyhedron_Copy(P), P->Dimension, options);
	if (!gf)
	    gf = P_gf;
	else {
	    gf->add_union(P_gf, options);
	    delete P_gf;
	}
    }
    /* we actually only need the convex union of the parameter space
     * but the reducer classes currently expect a polyhedron in
     * the combined space
     */
    Polyhedron_Free(gf->context);
    gf->context = DomainConvex(D2, options->MaxRays);

    gf2 = gf->summate(D2->Dimension - nparam, options);

    delete gf;
    if (D != D2)
	Domain_Free(D2);
    Domain_Free(D);
    return gf2;
}

gen_fun* barvinok_enumerate_union_series(Polyhedron *D, Polyhedron* C, 
					 unsigned MaxRays)
{
    gen_fun *gf;
    barvinok_options *options = barvinok_options_new_with_defaults();
    options->MaxRays = MaxRays;
    gf = barvinok_enumerate_union_series_with_options(D, C, options);
    barvinok_options_free(options);
    return gf;
}
