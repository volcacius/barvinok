#include <barvinok/barvinok.h>
#include <omega.h>
#include <vector>
#include "count.h"
#include "config.h"

#include <iostream>
using namespace std;

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

typedef std::vector<Variable_ID> varvector;

static void max_index(Constraint_Handle c, varvector& vv, varvector& params)
{
    for (Constr_Vars_Iter cvi(c); cvi; ++cvi) {
	Variable_ID var = (*cvi).var;
	if (find(vv.begin(), vv.end(), var) == vv.end() &&
	    find(params.begin(), params.end(), var) == params.end())
		vv.push_back(var);
    }
}

static void set_constraint(Matrix *M, int row,
			   Constraint_Handle c, varvector& vv, int geq)
{
    value_set_si(M->p[row][0], geq);
    for (int i = 0; i < vv.size(); ++i)
	value_set_si(M->p[row][1+i], c.get_coef(vv[i]));
    value_set_si(M->p[row][1+vv.size()], c.get_const());
}


evalue *count_relation(Relation& r)
{
    r.simplify();

    varvector vv;
    varvector params;
    if (r.is_set())
	for (int j = 1; j <= r.n_set(); ++j)
	    vv.push_back(r.set_var(j));
    else {
	for (int j = 1; j <= r.n_inp(); ++j)
	    vv.push_back(r.input_var(j));
	for (int j = 1; j <= r.n_out(); ++j)
	    vv.push_back(r.output_var(j));
    }
    int dim = vv.size();

    const Variable_ID_Tuple * globals = r.global_decls();
    for (int i = 0; i < globals->size(); ++i)
	params.push_back(r.get_local((*globals)[i+1]));

    int d = 0;
    evalue *EP = NULL;
    for (DNF_Iterator di(r.query_DNF()); di; ++di, ++d) {
	assert(d == 0);

	int c = 0;
	for (EQ_Iterator ei = (*di)->EQs(); ei; ++ei, ++c)
	    max_index((*ei), vv, params);
	for (GEQ_Iterator gi = (*di)->GEQs(); gi; ++gi, ++c)
	    max_index((*gi), vv, params);
	for (int i = 0; i < params.size(); ++i)
	    vv.push_back(params[i]);

	Matrix *M = Matrix_Alloc(c, vv.size() + 2);
	int row = 0;
	for (EQ_Iterator ei = (*di)->EQs(); ei; ++ei)
	    set_constraint(M, row++, (*ei), vv, 0);
	for (GEQ_Iterator gi = (*di)->GEQs(); gi; ++gi)
	    set_constraint(M, row++, (*gi), vv, 1);
	Polyhedron *P = Constraints2Polyhedron(M, MAXRAYS);
	Matrix_Free(M);
	int exist = P->Dimension - params.size() - dim;
	EP = barvinok_enumerate_e(P, exist, params.size(), MAXRAYS);
	Polyhedron_Free(P);
    }
    return EP;
}
