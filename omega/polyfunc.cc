#include <assert.h>
#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/options.h>
#include <barvinok/bernstein.h>
#include "ex_convert.h"
#include "polyfunc.h"
#include "param_util.h"
#include "omega/convert.h"

using namespace GiNaC;
using namespace bernstein;
using namespace std;

#define ALLOC(type) (type*)malloc(sizeof(type))

static void extract_vars(Map<Variable_Ref *, GiNaC::ex>& variableMap,
			 const varvector &vv, exvector &exvars)
{
    for (int i = 0; i < vv.size(); ++i) {
	ex var = 0;
	foreach_map(vr,Variable_Ref *,s,ex,variableMap, {
	    if (vr->vid == vv[i]) {
		var = s;
		break;
	    }
	});
	assert(var != 0);
	exvars.push_back(var);
    }
}

static void extract_params(Map<Variable_Ref *, GiNaC::ex>& variableMap,
			     const varvector &params, exvector &exparams)
{
    for (int i = 0; i < params.size(); ++i) {
	Global_Var_ID global = params[i]->get_global_var();
	ex var = 0;
	foreach_map(vr,Variable_Ref *,s,ex,variableMap, {
	    if (vr->g == global) {
		var = s;
		break;
	    }
	});
	if (var != 0)
	    exparams.push_back(var);
	else
	    exparams.push_back(symbol(params[i]->char_name()));
    }
}

void maximize(PolyFunc *polyfunc, Map<Variable_Ref *, GiNaC::ex>& variableMap)
{
    varvector vv;
    varvector params;
    Param_Polyhedron *PP;
    exvector exvars;
    exvector exparams;
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    cout << "maximize " << polyfunc->poly << " over ";
    polyfunc->domain.simplify();
    polyfunc->domain.print(stdout);

    Polyhedron *D = relation2Domain(polyfunc->domain, vv, params,
					options->MaxRays);
    assert(!D->next);
    assert(polyfunc->domain.is_set());
    int dim = polyfunc->domain.n_set();
    assert(D->Dimension == dim + params.size());
    Polyhedron *ctx = Universe_Polyhedron(params.size());

    extract_vars(variableMap, vv, exvars);
    extract_params(variableMap, params, exparams);

    PP = Polyhedron2Param_Polyhedron(D, ctx, options);
    piecewise_lst *pl = new piecewise_lst(exparams);
    for (Param_Domain *Q = PP->D; Q; Q = Q->next) {
	GiNaC::matrix VM = domainVertices(PP, Q, exparams);
	lst coeffs = bernsteinExpansion(VM, polyfunc->poly, exvars, exparams);
	pl->list.push_back(guarded_lst(Polyhedron_Copy(Q->Domain), coeffs));
    }
    cout << "coefficients: " << *pl << endl;
    pl->maximize();
    cout << "maximum: " << *pl << endl;
    delete pl;
    Param_Polyhedron_Free(PP);

    Polyhedron_Free(ctx);
    Domain_Free(D);
    barvinok_options_free(options);
}

evalue *summate(PolyFunc *polyfunc, Map<Variable_Ref *, GiNaC::ex>& variableMap)
{
    varvector vv;
    varvector params;
    exvector exvars;
    exvector exparams;
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    cout << "sum " << polyfunc->poly << " over ";
    polyfunc->domain.simplify();
    polyfunc->domain.print(stdout);

    Polyhedron *D = relation2Domain(polyfunc->domain, vv, params,
					options->MaxRays);
    assert(!D->next);
    assert(polyfunc->domain.is_set());
    int dim = polyfunc->domain.n_set();
    assert(D->Dimension == dim + params.size());

    extract_vars(variableMap, vv, exvars);
    extract_params(variableMap, params, exparams);
    for (int i = 0; i < exparams.size(); ++i)
	exvars.push_back(exparams[i]);

    evalue *e = ALLOC(evalue);
    value_init(e->d);
    value_set_si(e->d, 0);
    e->x.p = new_enode(::partition, 2, D->Dimension);
    EVALUE_SET_DOMAIN(e->x.p->arr[0], D);
    value_clear(e->x.p->arr[1].d);
    e->x.p->arr[1] = *ex2evalue(polyfunc->poly, exvars);

    evalue *sum = barvinok_summate(e, dim, options);
    evalue_free(e);

    barvinok_options_free(options);
    return sum;
}
