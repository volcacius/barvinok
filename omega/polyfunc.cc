#include <barvinok/bernstein.h>
#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include "polyfunc.h"
#include "omega/convert.h"

using namespace GiNaC;
using namespace bernstein;
using namespace std;

void maximize(PolyFunc *polyfunc, Map<Variable_Ref *, GiNaC::ex>& variableMap)
{
    varvector vv;
    varvector params;
    Param_Polyhedron *PP;
    exvector exvars;
    exvector exparams;

    cout << "maximize " << polyfunc->poly << " over ";
    polyfunc->domain.print(stdout);

    Polyhedron *D = relation2Domain(polyfunc->domain, vv, params);
    assert(!D->next);
    assert(polyfunc->domain.is_set());
    int dim = polyfunc->domain.n_set();
    assert(D->Dimension == dim + params.size());
    Polyhedron *ctx = Universe_Polyhedron(params.size());

    for (int i = 0; i < dim; ++i) {
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
    for (int i = dim; i < vv.size(); ++i) {
	Global_Var_ID global = vv[i]->get_global_var();
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
	    exparams.push_back(symbol(vv[i]->char_name()));
    }

    PP = Polyhedron2Param_Domain(D, ctx, 0);
    piecewise_lst_s *pl = new piecewise_lst_s(exparams);
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
}
