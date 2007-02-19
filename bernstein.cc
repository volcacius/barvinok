#include <vector>
#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>
#include <barvinok/options.h>

using namespace GiNaC;
using namespace bernstein;

using std::vector;
using std::cerr;
using std::endl;

namespace barvinok {

ex evalue2ex(evalue *e, const exvector& vars)
{
    if (value_notzero_p(e->d))
	return value2numeric(e->x.n)/value2numeric(e->d);
    if (e->x.p->type != polynomial)
	return fail();
    ex poly = 0;
    for (int i = e->x.p->size-1; i >= 0; --i) {
	poly *= vars[e->x.p->pos-1];
	ex t = evalue2ex(&e->x.p->arr[i], vars);
	if (is_exactly_a<fail>(t))
	    return t;
	poly += t;
    }
    return poly;
}

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 : 
	  p->type == flooring ? 1 : 0;
}

static ex evalue2ex_r(const evalue *e, const exvector& vars,
		      exvector& floorvar, vector<const evalue *>& floors)
{
    if (value_notzero_p(e->d))
	return value2numeric(e->x.n)/value2numeric(e->d);
    ex base_var = 0;
    ex poly = 0;

    switch (e->x.p->type) {
    case polynomial:
	base_var = vars[e->x.p->pos-1];
	break;
    case flooring:
	for (int i = 0; i < floors.size(); ++i) {
	    if (eequal(&e->x.p->arr[0], floors[i])) {
		base_var = floorvar[i];
		break;
	    }
	}
	if (base_var != 0)
	    break;
	char name[20];
	snprintf(name, sizeof(name), "f%d", floors.size());
	floorvar.push_back(base_var = symbol(name));
	floors.push_back(&e->x.p->arr[0]);
	break;
    default:
	return fail();
    }

    int offset = type_offset(e->x.p);
    for (int i = e->x.p->size-1; i >= offset; --i) {
	poly *= base_var;
	ex t = evalue2ex_r(&e->x.p->arr[i], vars, floorvar, floors);
	if (is_exactly_a<fail>(t))
	    return t;
	poly += t;
    }
    return poly;
}

/* For each t = floor(e/d), set up two constraints
 *
 *	 e - d t       >= 0
 *	-e + d t + d-1 >= 0
 *
 * e is assumed to be an affine expression.
 */
static Matrix *setup_floor_constraints(const vector<const evalue*> floors, int nvar)
{
    int extra = floors.size();
    if (!extra)
	return NULL;
    Matrix *M = Matrix_Alloc(2*extra, 1+extra+nvar+1);
    for (int i = 0; i < extra; ++i) {
	Value *d = &M->p[2*i][1+i];
	value_set_si(*d, 1);
	evalue_denom(floors[i], d);
	const evalue *e;
	for (e = floors[i]; value_zero_p(e->d); e = &e->x.p->arr[0]) {
	    assert(e->x.p->type == polynomial);
	    assert(e->x.p->size == 2);
	    evalue *c = &e->x.p->arr[1];
	    value_multiply(M->p[2*i][1+extra+e->x.p->pos-1], *d, c->x.n);
	    value_division(M->p[2*i][1+extra+e->x.p->pos-1],
			   M->p[2*i][1+extra+e->x.p->pos-1], c->d);
	}
	value_multiply(M->p[2*i][1+extra+nvar], *d, e->x.n);
	value_division(M->p[2*i][1+extra+nvar], M->p[2*i][1+extra+nvar], e->d);
	value_oppose(*d, *d);
	value_set_si(M->p[2*i][0], -1);
	Vector_Scale(M->p[2*i], M->p[2*i+1], M->p[2*i][0], 1+extra+nvar+1);
	value_set_si(M->p[2*i][0], 1);
	value_subtract(M->p[2*i+1][1+extra+nvar], M->p[2*i+1][1+extra+nvar], *d);
	value_decrement(M->p[2*i+1][1+extra+nvar], M->p[2*i+1][1+extra+nvar]);
    }
    return M;
}

ex evalue2ex(const evalue *e, const exvector& vars, exvector& floorvar, Matrix **C)
{
    vector<const evalue *> floors;
    ex poly = evalue2ex_r(e, vars, floorvar, floors);
    assert(C);
    Matrix *M = setup_floor_constraints(floors, vars.size());
    *C = M;
    return poly;
}

/* if the evalue is a relation, we use the relation to cut off the 
 * the edges of the domain
 */
static void evalue_extract_poly(evalue *e, int i, Polyhedron **D, evalue **poly,
				unsigned MaxRays)
{
    *D = EVALUE_DOMAIN(e->x.p->arr[2*i]);
    *poly = e = &e->x.p->arr[2*i+1];
    if (value_notzero_p(e->d))
	return;
    if (e->x.p->type != relation)
	return;
    if (e->x.p->size > 2)
	return;
    evalue *fr = &e->x.p->arr[0];
    assert(value_zero_p(fr->d));
    assert(fr->x.p->type == fractional);
    assert(fr->x.p->size == 3);
    Matrix *T = Matrix_Alloc(2, (*D)->Dimension+1);
    value_set_si(T->p[1][(*D)->Dimension], 1);

    /* convert argument of fractional to polylib */
    /* the argument is assumed to be linear */
    evalue *p = &fr->x.p->arr[0];
    evalue_denom(p, &T->p[1][(*D)->Dimension]);
    for (;value_zero_p(p->d); p = &p->x.p->arr[0]) {
	assert(p->x.p->type == polynomial);
	assert(p->x.p->size == 2);
	assert(value_notzero_p(p->x.p->arr[1].d));
	int pos = p->x.p->pos - 1;
	value_assign(T->p[0][pos], p->x.p->arr[1].x.n);
	value_multiply(T->p[0][pos], T->p[0][pos], T->p[1][(*D)->Dimension]);
	value_division(T->p[0][pos], T->p[0][pos], p->x.p->arr[1].d);
    }
    int pos = (*D)->Dimension;
    value_assign(T->p[0][pos], p->x.n);
    value_multiply(T->p[0][pos], T->p[0][pos], T->p[1][(*D)->Dimension]);
    value_division(T->p[0][pos], T->p[0][pos], p->d);

    Polyhedron *E = NULL;
    for (Polyhedron *P = *D; P; P = P->next) {
	Polyhedron *I = Polyhedron_Image(P, T, MaxRays);
	I = DomainConstraintSimplify(I, MaxRays);
	Polyhedron *R = Polyhedron_Preimage(I, T, MaxRays);
	Polyhedron_Free(I);
	Polyhedron *next = P->next;
	P->next = NULL;
	Polyhedron *S = DomainIntersection(P, R, MaxRays);
	Polyhedron_Free(R);
	P->next = next;
	if (emptyQ2(S))
	    Polyhedron_Free(S);
	else
	    E = DomainConcat(S, E);
    }
    Matrix_Free(T);

    *D = E;
    *poly = &e->x.p->arr[1];
}

piecewise_lst *evalue_bernstein_coefficients(piecewise_lst *pl_all, evalue *e, 
				      Polyhedron *ctx, const exvector& params)
{
    piecewise_lst *pl;
    barvinok_options *options = barvinok_options_new_with_defaults();
    pl = evalue_bernstein_coefficients(pl_all, e, ctx, params, options);
    barvinok_options_free(options);
    return pl;
}

piecewise_lst *evalue_bernstein_coefficients(piecewise_lst *pl_all, evalue *e, 
				      Polyhedron *ctx, const exvector& params,
				      barvinok_options *options)
{
    unsigned nparam = ctx->Dimension;
    if (EVALUE_IS_ZERO(*e))
	return pl_all;
    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);
    assert(e->x.p->size >= 2);
    unsigned nvars = EVALUE_DOMAIN(e->x.p->arr[0])->Dimension - nparam;
    unsigned PP_MaxRays = options->MaxRays;
    if (PP_MaxRays & POL_NO_DUAL)
	PP_MaxRays = 0;

    exvector vars = constructVariableVector(nvars, "v");
    exvector allvars = vars;
    allvars.insert(allvars.end(), params.begin(), params.end());

    for (int i = 0; i < e->x.p->size/2; ++i) {
	Param_Polyhedron *PP;
	Polyhedron *E;
	evalue *EP;
	Matrix *M;
	exvector floorvar;

	evalue_extract_poly(e, i, &E, &EP, options->MaxRays);
	ex poly = evalue2ex(EP, allvars, floorvar, &M);
	floorvar.insert(floorvar.end(), vars.begin(), vars.end());
	if (M) {
	    Polyhedron *AE = align_context(E, M->NbColumns-2, options->MaxRays);
	    if (E != EVALUE_DOMAIN(e->x.p->arr[2*i]))
		Domain_Free(E);
	    E = DomainAddConstraints(AE, M, options->MaxRays);
	    Matrix_Free(M);
	    Domain_Free(AE);
	}
	if (is_exactly_a<fail>(poly)) {
	    if (E != EVALUE_DOMAIN(e->x.p->arr[2*i]))
		Domain_Free(E);
	    delete pl_all;
	    return NULL;
	}
	for (Polyhedron *P = E; P; P = P->next) {
	    Polyhedron *next = P->next;
	    piecewise_lst *pl = new piecewise_lst(params);
	    Polyhedron *P1 = P;
	    P->next = NULL;
	    PP = Polyhedron2Param_Domain(P, ctx, PP_MaxRays);
	    for (Param_Domain *Q = PP->D; Q; Q = Q->next) {
		matrix VM = domainVertices(PP, Q, params);
		lst coeffs = bernsteinExpansion(VM, poly, floorvar, params);
		pl->list.push_back(guarded_lst(Polyhedron_Copy(Q->Domain), coeffs));
	    }
	    Param_Polyhedron_Free(PP);
	    if (!pl_all)
		pl_all = pl;
	    else {
		pl_all->combine(*pl);
		delete pl;
	    }
	    P->next = next;
	}
	if (E != EVALUE_DOMAIN(e->x.p->arr[2*i]))
	    Domain_Free(E);
    }
    return pl_all;
}

}
