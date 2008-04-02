#include <assert.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include <omega.h>
#include "omega/convert.h"
#include "count.h"

#include <iostream>
using namespace std;

evalue *count_relation(Relation& r)
{
    varvector vv;
    varvector params;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    Polyhedron *D = relation2Domain(r, vv, params, options->MaxRays);
    int dim = r.is_set() ? r.n_set() : r.n_inp() + r.n_out();

    int d = 0;
    evalue *EP = NULL;
    for (Polyhedron *P = D; P; P = P->next, ++d) {
	assert(d == 0);
	int exist = P->Dimension - params.size() - dim;
	EP = barvinok_enumerate_e(P, exist, params.size(), options->MaxRays);
    }
    reduce_evalue(EP);
    Domain_Free(D);
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
