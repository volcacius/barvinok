#include <vector>
#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>
#include <barvinok/options.h>

using namespace GiNaC;
using namespace bernstein;

using std::pair;
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

typedef pair<bool, const evalue *> typed_evalue;

static ex evalue2ex_add_var(evalue *e, exvector& extravar,
			    vector<typed_evalue>& expr, bool is_fract)
{
    ex base_var = 0;

    for (int i = 0; i < expr.size(); ++i) {
	if (is_fract == expr[i].first && eequal(e, expr[i].second)) {
	    base_var = extravar[i];
	    break;
	}
    }
    if (base_var != 0)
	return base_var;

    char name[20];
    snprintf(name, sizeof(name), "f%c%d", is_fract ? 'r' : 'l', expr.size());
    extravar.push_back(base_var = symbol(name));
    expr.push_back(typed_evalue(is_fract, e));

    return base_var;
}

static ex evalue2ex_r(const evalue *e, const exvector& vars,
		      exvector& extravar, vector<typed_evalue>& expr,
		      Vector *coset)
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
	base_var = evalue2ex_add_var(&e->x.p->arr[0], extravar, expr, false);
	break;
    case fractional:
	base_var = evalue2ex_add_var(&e->x.p->arr[0], extravar, expr, true);
	break;
    case periodic:
	assert(coset);
	return evalue2ex_r(&e->x.p->arr[VALUE_TO_INT(coset->p[e->x.p->pos-1])],
			   vars, extravar, expr, coset);
    default:
	return fail();
    }

    int offset = type_offset(e->x.p);
    for (int i = e->x.p->size-1; i >= offset; --i) {
	poly *= base_var;
	ex t = evalue2ex_r(&e->x.p->arr[i], vars, extravar, expr, coset);
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
 *
 * For each t = fract(e/d), set up two constraints
 *
 *          -d t + d-1 >= 0
 *             t       >= 0
 */
static Matrix *setup_constraints(const vector<typed_evalue> expr, int nvar)
{
    int extra = expr.size();
    if (!extra)
	return NULL;
    Matrix *M = Matrix_Alloc(2*extra, 1+extra+nvar+1);
    for (int i = 0; i < extra; ++i) {
	Value *d = &M->p[2*i][1+i];
	value_set_si(*d, 1);
	evalue_denom(expr[i].second, d);
	if (expr[i].first) {
	    value_set_si(M->p[2*i][0], 1);
	    value_decrement(M->p[2*i][1+extra+nvar], *d);
	    value_oppose(*d, *d);
	    value_set_si(M->p[2*i+1][0], 1);
	    value_set_si(M->p[2*i+1][1+i], 1);
	} else {
	    const evalue *e;
	    for (e = expr[i].second; value_zero_p(e->d); e = &e->x.p->arr[0]) {
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
    }
    return M;
}

static bool evalue_is_periodic(const evalue *e, Vector *periods)
{
    int i, offset;
    bool is_periodic = false;

    if (value_notzero_p(e->d))
	return false;

    assert(e->x.p->type != partition);
    if (e->x.p->type == periodic) {
	Value size;
	value_init(size);
	value_set_si(size, e->x.p->size);
	value_lcm(periods->p[e->x.p->pos-1], size, &periods->p[e->x.p->pos-1]);
	value_clear(size);
	is_periodic = true;
    }
    offset = type_offset(e->x.p);
    for (i = e->x.p->size-1; i >= offset; --i)
	is_periodic = evalue_is_periodic(&e->x.p->arr[i], periods) || is_periodic;
    return is_periodic;
}

static ex evalue2lst(const evalue *e, const exvector& vars,
		     exvector& extravar, vector<typed_evalue>& expr,
		     Vector *periods)
{
    Vector *coset = Vector_Alloc(periods->Size);
    lst list;
    for (;;) {
	int i;
	list.append(evalue2ex_r(e, vars, extravar, expr, coset));
	for (i = coset->Size-1; i >= 0; --i) {
	    value_increment(coset->p[i], coset->p[i]);
	    if (value_lt(coset->p[i], periods->p[i]))
		break;
	    value_set_si(coset->p[i], 0);
	}
	if (i < 0)
	    break;
    }
    Vector_Free(coset);
    return list;
}

ex evalue2ex(const evalue *e, const exvector& vars, exvector& floorvar,
	     Matrix **C, Vector **p)
{
    vector<typed_evalue> expr;
    Vector *periods = Vector_Alloc(vars.size());
    assert(p);
    assert(C);
    for (int i = 0; i < periods->Size; ++i)
	value_set_si(periods->p[i], 1);
    if (evalue_is_periodic(e, periods)) {
	*p = periods;
	*C = NULL;
	lst list;
	return list;
    } else {
	Vector_Free(periods);
	*p = NULL;
	ex poly = evalue2ex_r(e, vars, floorvar, expr, NULL);
	Matrix *M = setup_constraints(expr, vars.size());
	*C = M;
	return poly;
    }
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

static piecewise_lst *bernstein_coefficients(piecewise_lst *pl_all,
			    Polyhedron *D, const ex& poly,
			    Polyhedron *ctx,
			    const exvector& params, const exvector& floorvar,
			    barvinok_options *options)
{
    if (emptyQ2(D))
	return pl_all;

    unsigned PP_MaxRays = options->MaxRays;
    if (PP_MaxRays & POL_NO_DUAL)
	PP_MaxRays = 0;

    for (Polyhedron *P = D; P; P = P->next) {
	Param_Polyhedron *PP;
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
    return pl_all;
}

/* Compute the coefficients of the polynomial corresponding to each coset
 * on its own domain.  This allows us to cut the domain on multiples of
 * the period.
 * To perform the cutting for a coset "i mod n = c" we map the domain
 * to the quotient space trough "i = i' n + c", simplify the constraints
 * (implicitly) and then map back to the original space.
 */
static piecewise_lst *bernstein_coefficients_periodic(piecewise_lst *pl_all,
			    Polyhedron *D, const evalue *e,
			    Polyhedron *ctx, const exvector& vars,
			    const exvector& params, Vector *periods,
			    barvinok_options *options)
{
    assert(D->Dimension == periods->Size);
    Matrix *T = Matrix_Alloc(D->Dimension+1, D->Dimension+1);
    Matrix *T2 = Matrix_Alloc(D->Dimension+1, D->Dimension+1);
    Vector *coset = Vector_Alloc(periods->Size);
    exvector extravar;
    vector<typed_evalue> expr;
    exvector allvars = vars;
    allvars.insert(allvars.end(), params.begin(), params.end());

    value_set_si(T2->p[D->Dimension][D->Dimension], 1);
    for (int i = 0; i < D->Dimension; ++i) {
	value_assign(T->p[i][i], periods->p[i]);
	value_lcm(T2->p[D->Dimension][D->Dimension], periods->p[i],
		  &T2->p[D->Dimension][D->Dimension]);
    }
    value_set_si(T->p[D->Dimension][D->Dimension], 1);
    for (int i = 0; i < D->Dimension; ++i)
	value_division(T2->p[i][i], T2->p[D->Dimension][D->Dimension],
		       periods->p[i]);
    for (;;) {
	int i;
	ex poly = evalue2ex_r(e, allvars, extravar, expr, coset);
	assert(extravar.size() == 0);
	assert(expr.size() == 0);
	Polyhedron *E = DomainPreimage(D, T, options->MaxRays);
	Polyhedron *F = DomainPreimage(E, T2, options->MaxRays);
	Polyhedron_Free(E);
	if (!emptyQ2(F))
	    pl_all = bernstein_coefficients(pl_all, F, poly, ctx, params,
					    vars, options);
	Polyhedron_Free(F);
	for (i = D->Dimension-1; i >= 0; --i) {
	    value_increment(coset->p[i], coset->p[i]);
	    value_increment(T->p[i][D->Dimension], T->p[i][D->Dimension]);
	    value_subtract(T2->p[i][D->Dimension], T2->p[i][D->Dimension],
			   T2->p[i][i]);
	    if (value_lt(coset->p[i], periods->p[i]))
		break;
	    value_set_si(coset->p[i], 0);
	    value_set_si(T->p[i][D->Dimension], 0);
	    value_set_si(T2->p[i][D->Dimension], 0);
	}
	if (i < 0)
	    break;
    }
    Vector_Free(coset);
    Matrix_Free(T);
    Matrix_Free(T2);
    return pl_all;
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

    exvector vars = constructVariableVector(nvars, "v");
    exvector allvars = vars;
    allvars.insert(allvars.end(), params.begin(), params.end());

    for (int i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *E;
	evalue *EP;
	Matrix *M;
	Vector *periods;
	exvector floorvar;

	evalue_extract_poly(e, i, &E, &EP, options->MaxRays);
	ex poly = evalue2ex(EP, allvars, floorvar, &M, &periods);
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
	if (periods)
	    pl_all = bernstein_coefficients_periodic(pl_all, E, EP, ctx, vars,
						     params, periods, options);
	else
	    pl_all = bernstein_coefficients(pl_all, E, poly, ctx, params,
					    floorvar, options);
	if (periods)
	    Vector_Free(periods);
	if (E != EVALUE_DOMAIN(e->x.p->arr[2*i]))
	    Domain_Free(E);
    }
    return pl_all;
}

}
