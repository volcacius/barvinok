#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

using namespace GiNaC;
using namespace bernstein;

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

/* if the evalue is a relation, we use the relation to cut off the 
 * the edges of the domain
 */
static void evalue_extract_poly(evalue *e, int i, Polyhedron **D, evalue **poly)
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
	Polyhedron *I = Polyhedron_Image(P, T, MAXRAYS);
	I = DomainConstraintSimplify(I, MAXRAYS);
	Polyhedron *R = Polyhedron_Preimage(I, T, MAXRAYS);
	Polyhedron_Free(I);
	Polyhedron *next = P->next;
	P->next = NULL;
	Polyhedron *S = DomainIntersection(P, R, MAXRAYS);
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
    unsigned nparam = ctx->Dimension;
    if (EVALUE_IS_ZERO(*e))
	return pl_all;
    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);
    assert(e->x.p->size >= 2);
    unsigned nvars = EVALUE_DOMAIN(e->x.p->arr[0])->Dimension - nparam;

    exvector vars = constructVariableVector(nvars, "v");
    exvector allvars = vars;
    allvars.insert(allvars.end(), params.begin(), params.end());

    for (int i = 0; i < e->x.p->size/2; ++i) {
	Param_Polyhedron *PP;
	Polyhedron *E;
	evalue *EP;
	evalue_extract_poly(e, i, &E, &EP);
	ex poly = evalue2ex(EP, allvars);
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
	    PP = Polyhedron2Param_Domain(P, ctx, 0);
	    for (Param_Domain *Q = PP->D; Q; Q = Q->next) {
		matrix VM = domainVertices(PP, Q, params);
		lst coeffs = bernsteinExpansion(VM, poly, vars, params);
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
