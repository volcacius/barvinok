#include <vector>
#include <assert.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include <omega.h>
#include "omega/convert.h"
#include "normalization.h"
#include "count.h"
#include "config.h"

#include <iostream>
using namespace std;

#define ALLOC(t,p) p = (t*)malloc(sizeof(*p))

#ifdef USE_PARKER

#include "parker/count_solutions.h"

/*
 * Use parker's method to compute the number of integer points in D.
 * Since this method assumes all variables are non-negative,
 * we have to transform the input polytope first.
 */
evalue *barvinok_enumerate_parker(Polyhedron *D,
					unsigned nvar, unsigned nparam,
					unsigned MaxRays)
{
    Polyhedron *R;
    evalue *res;

    if (nparam != 0) {
	fprintf(stderr, "parker method doesn't support parameters\n");
	return NULL;
    }
    R = skew_to_positive_orthant(D, nvar, MaxRays);
    Relation r = Domain2relation(R, nvar, 0, NULL);
    Polyhedron_Free(R);
    double d = count_solutions(r);
    ALLOC(evalue, res);
    value_init(res->d);
    value_set_si(res->d, 0);
    res->x.p = new_enode(::partition, 2, 0);
    EVALUE_SET_DOMAIN(res->x.p->arr[0], Universe_Polyhedron(0));
    value_set_si(res->x.p->arr[1].d, 1);
    value_init(res->x.p->arr[1].x.n);
    value_set_double(res->x.p->arr[1].x.n, d);
    return res;
}

#else

evalue *barvinok_enumerate_parker(Polyhedron *D,
					unsigned nvar, unsigned nparam,
					unsigned MaxRays)
{
    fprintf(stderr, "support for parker method not compiled in\n");
    return NULL;
}

#endif

static evalue *count_relation_barvinok(Polyhedron *D, unsigned nvar,
					unsigned nparam, unsigned MaxRays)
{
    int d = 0;
    evalue *EP = NULL;
    for (Polyhedron *P = D; P; P = P->next, ++d) {
	assert(d == 0);
	int exist = P->Dimension - nparam - nvar;
	EP = barvinok_enumerate_e(P, exist, nparam, MaxRays);
    }
    reduce_evalue(EP);
    Domain_Free(D);
    return EP;
}

evalue *count_relation(Relation& r, int counting_method)
{
    varvector vv;
    varvector params;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    Polyhedron *D = relation2Domain(r, vv, params, options->MaxRays);

    evalue *EP = NULL;
    if (counting_method == COUNT_RELATION_BARVINOK)
	EP = count_relation_barvinok(D, vv.size(), params.size(),
					options->MaxRays);
    else if (counting_method == COUNT_RELATION_PARKER)
	EP = barvinok_enumerate_parker(D, vv.size(), params.size(),
					options->MaxRays);
    barvinok_options_free(options);
    return EP;
}

evalue *rank_relation(Relation& r)
{
    varvector vv;
    varvector params;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    Polyhedron *D = relation2Domain(r, vv, params, options->MaxRays);
    int dim = r.is_set() ? r.n_set() : r.n_inp() + r.n_out();

    evalue *EP = NULL;
    if (D) {
	assert(D->next == NULL);
	Polyhedron *C = Universe_Polyhedron(params.size());
	EP = barvinok_lexsmaller_ev(D, D, dim, C, options->MaxRays);
	Polyhedron_Free(C);
    }
    Domain_Free(D);
    barvinok_options_free(options);
    return EP;
}

evalue *count_lexsmaller(Relation& r, Relation& domain)
{
    varvector P_vv;
    varvector P_params;
    varvector D_vv;
    varvector D_params;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    Polyhedron *P = relation2Domain(r, P_vv, P_params, options->MaxRays);
    int P_dim = r.is_set() ? r.n_set() : r.n_inp() + r.n_out();
    Polyhedron *D = relation2Domain(domain, D_vv, D_params, options->MaxRays);
    int D_dim = r.is_set() ? r.n_set() : r.n_inp() + r.n_out();
    assert(P_dim == D_dim);

    evalue *EP = NULL;
    if (P && D) {
	assert(P->next == NULL);
	assert(D->next == NULL);
	Polyhedron *C = Universe_Polyhedron(D_params.size());
	EP = barvinok_lexsmaller_ev(P, D, D_dim, C, options->MaxRays);
	Polyhedron_Free(C);
    }
    Domain_Free(P);
    Domain_Free(D);
    barvinok_options_free(options);
    return EP;
}
