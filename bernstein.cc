#include <assert.h>
#include <vector>
#include <bernstein/bernstein.h>
#include <bernstein/piecewise_lst.h>
#include <isl_set_polylib.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>
#include <barvinok/options.h>
#include "reduce_domain.h"

using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

using std::pair;
using std::vector;
using std::cerr;
using std::endl;

namespace barvinok {

ex evalue2ex(evalue *e, const exvector& vars)
{
    if (value_pos_p(e->d))
	return value2numeric(e->x.n)/value2numeric(e->d);
    if (EVALUE_IS_NAN(*e))
	return fail();
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

/* For the argument e=(f/d) of a fractional, return (d-1)/d times
 * a variable in [0,1] (see setup_constraints).
 */
static ex evalue2ex_get_fract(evalue *e, exvector& extravar,
				vector<typed_evalue>& expr)
{
    ex f;
    Value d;
    ex den;
    value_init(d);
    value_set_si(d, 1);
    evalue_denom(e, &d);
    den = value2numeric(d);
    value_clear(d);
    f = (den-1)/den;

    ex base_var = evalue2ex_add_var(e, extravar, expr, true);
    base_var *= f;
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
    int rem;

    switch (e->x.p->type) {
    case polynomial:
	base_var = vars[e->x.p->pos-1];
	break;
    case flooring:
	base_var = evalue2ex_add_var(&e->x.p->arr[0], extravar, expr, false);
	break;
    case fractional:
	base_var = evalue2ex_get_fract(&e->x.p->arr[0], extravar, expr);
	break;
    case periodic:
	assert(coset);
	rem = VALUE_TO_INT(coset->p[e->x.p->pos-1]) % e->x.p->size;
	return evalue2ex_r(&e->x.p->arr[rem], vars, extravar, expr, coset);
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
	if (expr[i].first) {
	    value_set_si(M->p[2*i][0], 1);
	    value_set_si(M->p[2*i][1+i], -1);
	    value_set_si(M->p[2*i][1+extra+nvar], 1);
	    value_set_si(M->p[2*i+1][0], 1);
	    value_set_si(M->p[2*i+1][1+i], 1);
	} else {
	    Value *d = &M->p[2*i][1+i];
	    evalue_extract_affine(expr[i].second, M->p[2*i]+1+extra,
				  M->p[2*i]+1+extra+nvar, d);
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
	value_lcm(periods->p[e->x.p->pos-1], periods->p[e->x.p->pos-1], size);
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
static Polyhedron *relation_domain(Polyhedron *D, evalue *fr, unsigned MaxRays)
{
    assert(value_zero_p(fr->d));
    assert(fr->x.p->type == fractional);
    assert(fr->x.p->size == 3);
    Matrix *T = Matrix_Alloc(2, D->Dimension+1);
    value_set_si(T->p[1][D->Dimension], 1);

    /* convert argument of fractional to polylib */
    /* the argument is assumed to be linear */
    evalue *p = &fr->x.p->arr[0];
    evalue_denom(p, &T->p[1][D->Dimension]);
    for (;value_zero_p(p->d); p = &p->x.p->arr[0]) {
	assert(p->x.p->type == polynomial);
	assert(p->x.p->size == 2);
	assert(value_notzero_p(p->x.p->arr[1].d));
	int pos = p->x.p->pos - 1;
	value_assign(T->p[0][pos], p->x.p->arr[1].x.n);
	value_multiply(T->p[0][pos], T->p[0][pos], T->p[1][D->Dimension]);
	value_division(T->p[0][pos], T->p[0][pos], p->x.p->arr[1].d);
    }
    int pos = D->Dimension;
    value_assign(T->p[0][pos], p->x.n);
    value_multiply(T->p[0][pos], T->p[0][pos], T->p[1][D->Dimension]);
    value_division(T->p[0][pos], T->p[0][pos], p->d);

    Polyhedron *E = NULL;
    for (Polyhedron *P = D; P; P = P->next) {
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

    return E;
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
			    barvinok_options *options);

/* Recursively apply Bernstein expansion on P, optimizing over dims[i]
 * variables in each level.  The context ctx is assumed to have been adapted
 * to the first level in the recursion.
 */
static piecewise_lst *bernstein_coefficients_recursive(piecewise_lst *pl_all,
			    Polyhedron *P, const vector<int>& dims, const ex& poly,
			    Polyhedron *ctx,
			    const exvector& params, const exvector& vars,
			    barvinok_options *options)
{
    assert(dims.size() > 0);
    assert(ctx->Dimension == P->Dimension - dims[0]);
    piecewise_lst *pl;
    unsigned done = 0;
    for (int j = 0; j < dims.size(); ++j) {
	exvector pl_vars;
	pl_vars.insert(pl_vars.end(), vars.begin()+done, vars.begin()+done+dims[j]);
	exvector pl_params;
	pl_params.insert(pl_params.end(), vars.begin()+done+dims[j], vars.end());
	pl_params.insert(pl_params.end(), params.begin(), params.end());

	if (!j)
	    pl = bernstein_coefficients(NULL, P, poly, ctx,
					pl_params, pl_vars, options);
	else {
	    piecewise_lst *new_pl = NULL;
	    Polyhedron *U = Universe_Polyhedron(pl_params.size());

	    for (int i = 0; i < pl->list.size(); ++i) {
		Polyhedron *D = pl->list[i].first;
		lst polys = pl->list[i].second;
		new_pl = bernstein_coefficients(new_pl, D, polys, U, pl_params,
						pl_vars, options);
	    }

	    Polyhedron_Free(U);

	    delete pl;
	    pl = new_pl;
	}

	done += dims[j];

	if (!pl)
	    return pl_all;
    }

    if (!pl_all)
	pl_all = pl;
    else {
	pl_all->combine(*pl);
	delete pl;
    }

    return pl_all;
}

static piecewise_lst *bernstein_coefficients_full_recurse(piecewise_lst *pl_all,
			    Polyhedron *P, const ex& poly,
			    Polyhedron *ctx,
			    const exvector& params, const exvector& vars,
			    barvinok_options *options)
{
    Polyhedron *CR = align_context(ctx, P->Dimension-1, options->MaxRays);
    vector<int> dims(vars.size());
    for (int i = 0; i < dims.size(); ++i)
	dims[i] = 1;
    pl_all = bernstein_coefficients_recursive(pl_all, P, dims, poly, CR,
					      params, vars, options);
    Polyhedron_Free(CR);

    return pl_all;
}

static piecewise_lst *bernstein_coefficients_product(piecewise_lst *pl_all,
			    Polyhedron *F, Matrix *T, const ex& poly,
			    Polyhedron *ctx,
			    const exvector& params, const exvector& vars,
			    barvinok_options *options)
{
    if (emptyQ2(ctx))
	return pl_all;
    for (Polyhedron *G = F; G; G = G->next)
	if (emptyQ2(G))
	    return pl_all;

    unsigned nparam = params.size();
    unsigned nvar = vars.size();
    unsigned constraints;
    unsigned factors;
    Polyhedron *C = NULL;

    /* More context constraints */
    if (F->Dimension == ctx->Dimension) {
	C = F;
	F = F->next;
    }
    assert(F);
    assert(F->next);

    Matrix *M;
    Polyhedron *P;
    Polyhedron *PC;
    M = Matrix_Alloc(F->NbConstraints, 1+nvar+nparam+1);
    for (int i = 0; i < F->NbConstraints; ++i) {
	Vector_Copy(F->Constraint[i], M->p[i], 1+F->Dimension-nparam);
	Vector_Copy(F->Constraint[i]+1+F->Dimension-nparam,
		    M->p[i]+1+nvar, nparam+1);
    }
    P = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);

    factors = 1;
    constraints = C ? C->NbConstraints : 0;
    constraints += ctx->NbConstraints;
    for (Polyhedron *G = F->next; G; G = G->next) {
	constraints += G->NbConstraints;
	++factors;
    }

    unsigned total_var = nvar-(F->Dimension-nparam);
    unsigned skip = 0;
    unsigned c = 0;
    M = Matrix_Alloc(constraints, 1+total_var+nparam+1);
    for (Polyhedron *G = F->next; G; G = G->next) {
	unsigned this_var = G->Dimension - nparam;
	for (int i = 0; i < G->NbConstraints; ++i) {
	    value_assign(M->p[c+i][0], G->Constraint[i][0]);
	    Vector_Copy(G->Constraint[i]+1, M->p[c+i]+1+skip, this_var);
	    Vector_Copy(G->Constraint[i]+1+this_var, M->p[c+i]+1+total_var,
			nparam+1);
	}
	c += G->NbConstraints;
	skip += this_var;
    }
    assert(skip == total_var);
    if (C) {
	for (int i = 0; i < C->NbConstraints; ++i) {
	    value_assign(M->p[c+i][0], C->Constraint[i][0]);
	    Vector_Copy(C->Constraint[i]+1, M->p[c+i]+1+total_var,
			nparam+1);
	}
	c += C->NbConstraints;
    }
    for (int i = 0; i < ctx->NbConstraints; ++i) {
	value_assign(M->p[c+i][0], ctx->Constraint[i][0]);
	Vector_Copy(ctx->Constraint[i]+1, M->p[c+i]+1+total_var, nparam+1);
    }
    PC = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);

    exvector newvars = constructVariableVector(nvar, "t");
    matrix subs(1, nvar);
    for (int i = 0; i < nvar; ++i)
	for (int j = 0; j < nvar; ++j)
	    subs(0,i) += value2numeric(T->p[i][j]) * newvars[j];

    ex newpoly = replaceVariablesInPolynomial(poly, vars, subs);

    vector<int> dims(factors);
    for (int i = 0; F; ++i, F = F->next)
	dims[i] = F->Dimension-nparam;

    pl_all = bernstein_coefficients_recursive(pl_all, P, dims, newpoly, PC,
					      params, newvars, options);

    Polyhedron_Free(P);
    Polyhedron_Free(PC);

    return pl_all;
}

static piecewise_lst *bernstein_coefficients_polyhedron(piecewise_lst *pl_all,
			    Polyhedron *P, const ex& poly,
			    Polyhedron *ctx,
			    const exvector& params, const exvector& floorvar,
			    barvinok_options *options)
{
    if (Polyhedron_is_unbounded(P, ctx->Dimension, options->MaxRays)) {
	fprintf(stderr, "warning: unbounded domain skipped\n");
	Polyhedron_Print(stderr, P_VALUE_FMT, P);
	return pl_all;
    }

    if (options->bernstein_recurse & BV_BERNSTEIN_FACTORS) {
	Matrix *T = NULL;
	Polyhedron *F = Polyhedron_Factor(P, ctx->Dimension, &T, options->MaxRays);
	if (F) {
	    pl_all = bernstein_coefficients_product(pl_all, F, T, poly, ctx, params,
						    floorvar, options);
	    Domain_Free(F);
	    Matrix_Free(T);
	    return pl_all;
	}
    }
    if (floorvar.size() > 1 &&
		options->bernstein_recurse & BV_BERNSTEIN_INTERVALS)
	return bernstein_coefficients_full_recurse(pl_all, P, poly, ctx,
						   params, floorvar, options);

    unsigned PP_MaxRays = options->MaxRays;
    if (PP_MaxRays & POL_NO_DUAL)
	PP_MaxRays = 0;

    Param_Polyhedron *PP = Polyhedron2Param_Domain(P, ctx, PP_MaxRays);
    if (!PP)
	return pl_all;
    piecewise_lst *pl = new piecewise_lst(params, options->bernstein_optimize);

    int nd = -1;
    Polyhedron *TC = true_context(P, ctx, options->MaxRays);
    FORALL_REDUCED_DOMAIN(PP, TC, nd, options, i, PD, rVD)
	matrix VM = domainVertices(PP, PD, params);
	lst coeffs = bernsteinExpansion(VM, poly, floorvar, params);
	pl->add_guarded_lst(rVD, coeffs);
    END_FORALL_REDUCED_DOMAIN
    Polyhedron_Free(TC);

    Param_Polyhedron_Free(PP);
    if (!pl_all)
	pl_all = pl;
    else {
	pl_all->combine(*pl);
	delete pl;
    }

    return pl_all;
}

static piecewise_lst *bernstein_coefficients(piecewise_lst *pl_all,
			    Polyhedron *D, const ex& poly,
			    Polyhedron *ctx,
			    const exvector& params, const exvector& floorvar,
			    barvinok_options *options)
{
    if (!D->next && emptyQ2(D))
	return pl_all;

    for (Polyhedron *P = D; P; P = P->next) {
	/* This shouldn't happen */
	if (emptyQ2(P))
	    continue;
	Polyhedron *next = P->next;
	P->next = NULL;
	pl_all = bernstein_coefficients_polyhedron(pl_all, P, poly, ctx,
						   params, floorvar, options);
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
	value_lcm(T2->p[D->Dimension][D->Dimension],
		  T2->p[D->Dimension][D->Dimension], periods->p[i]);
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
	if (!emptyQ(F))
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

piecewise_lst *bernstein_coefficients_relation(piecewise_lst *pl_all,
		    Polyhedron *D, evalue *EP, Polyhedron *ctx,
		    const exvector& allvars, const exvector& vars,
		    const exvector& params, barvinok_options *options)
{
    if (value_zero_p(EP->d) && EP->x.p->type == relation) {
	Polyhedron *E = relation_domain(D, &EP->x.p->arr[0], options->MaxRays);
	if (E) {
	    pl_all = bernstein_coefficients_relation(pl_all, E, &EP->x.p->arr[1],
						     ctx, allvars, vars, params,
						     options);
	    Domain_Free(E);
	}
	/* In principle, we could cut off the edges of this domain too */
	if (EP->x.p->size > 2)
	    pl_all = bernstein_coefficients_relation(pl_all, D, &EP->x.p->arr[2],
						     ctx, allvars, vars, params,
						     options);
	return pl_all;
    }

    Matrix *M;
    exvector floorvar;
    Vector *periods;
    ex poly = evalue2ex(EP, allvars, floorvar, &M, &periods);
    floorvar.insert(floorvar.end(), vars.begin(), vars.end());
    Polyhedron *E = D;
    if (M) {
	Polyhedron *AE = align_context(D, M->NbColumns-2, options->MaxRays);
	E = DomainAddConstraints(AE, M, options->MaxRays);
	Matrix_Free(M);
	Domain_Free(AE);
    }
    if (is_exactly_a<fail>(poly)) {
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
    if (D != E)
	Domain_Free(E);

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
	pl_all = bernstein_coefficients_relation(pl_all,
			EVALUE_DOMAIN(e->x.p->arr[2*i]), &e->x.p->arr[2*i+1],
			ctx, allvars, vars, params, options);
    }
    return pl_all;
}

static __isl_give isl_qpolynomial *qp_from_ex(__isl_take isl_dim *dim,
	const GiNaC::ex ex, const GiNaC::exvector &params, int i)
{
	isl_qpolynomial *qp;
	isl_qpolynomial *base;
	int deg;
	int j;

	if (is_a<fail>(ex))
		return isl_qpolynomial_nan(dim);

	if (is_a<numeric>(ex)) {
		numeric r = ex_to<numeric>(ex);
		isl_int n;
		isl_int d;
		isl_int_init(n);
		isl_int_init(d);
		numeric2value(r.numer(), n);
		numeric2value(r.denom(), d);
		qp = isl_qpolynomial_rat_cst(dim, n, d);
		isl_int_clear(n);
		isl_int_clear(d);
		return qp;
	}

	deg = ex.degree(params[i]);
	if (deg == 0)
		return qp_from_ex(dim, ex, params, i + 1);

	base = isl_qpolynomial_var(isl_dim_copy(dim), isl_dim_param, i);
	qp = qp_from_ex(isl_dim_copy(dim), ex.coeff(params[i], deg),
			params, i + 1);

	for (j = deg - 1; j >= 0; --j) {
		qp = isl_qpolynomial_mul(qp, isl_qpolynomial_copy(base));
		qp = isl_qpolynomial_add(qp,
			    qp_from_ex(isl_dim_copy(dim), ex.coeff(params[i], j),
					params, i + 1));
	}

	isl_qpolynomial_free(base);
	isl_dim_free(dim);

	return qp;
}

__isl_give isl_qpolynomial *isl_qpolynomial_from_ginac(__isl_take isl_dim *dim,
	const GiNaC::ex &ex, const GiNaC::exvector &params)
{
	GiNaC::ex exp;

	if (!dim)
		return NULL;

	exp = ex.expand();
	return qp_from_ex(dim, exp, params, 0);
error:
	isl_dim_free(dim);
	return NULL;
}

__isl_give isl_qpolynomial_fold *isl_qpolynomial_fold_from_ginac(
	__isl_take isl_dim *dim, enum isl_fold type, const GiNaC::lst &lst,
	const GiNaC::exvector &params)
{
	isl_qpolynomial_fold *fold;
	lst::const_iterator j;

	fold = isl_qpolynomial_fold_empty(type, isl_dim_copy(dim));
	for (j = lst.begin(); j != lst.end(); ++j) {
		isl_qpolynomial *qp;
		isl_qpolynomial_fold *fold_i;

		qp = isl_qpolynomial_from_ginac(isl_dim_copy(dim), *j, params);
		fold_i = isl_qpolynomial_fold_alloc(type, qp);
		fold = isl_qpolynomial_fold_fold(fold, fold_i);
	}
	isl_dim_free(dim);
	return fold;
}

__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_fold_from_ginac(
	__isl_take isl_dim *dim, bernstein::piecewise_lst *pl,
	const GiNaC::exvector &params)
{
	enum isl_fold type;
	isl_pw_qpolynomial_fold *pwf;

	pwf = isl_pw_qpolynomial_fold_zero(isl_dim_copy(dim));

	type = pl->sign > 0 ? isl_fold_max
			    : pl->sign < 0 ? isl_fold_min : isl_fold_list;

	for (int i = 0; i < pl->list.size(); ++i) {
		isl_pw_qpolynomial_fold *pwf_i;
		isl_qpolynomial_fold *fold;
		isl_set *set;

		set = isl_set_new_from_polylib(pl->list[i].first,
				isl_dim_copy(dim));
		fold = isl_qpolynomial_fold_from_ginac(isl_dim_copy(dim),
				type, pl->list[i].second, params);
		pwf_i = isl_pw_qpolynomial_fold_alloc(set, fold);
		pwf = isl_pw_qpolynomial_fold_add_disjoint(pwf, pwf_i);
	}

	isl_dim_free(dim);

	return pwf;
}

}

struct isl_bound {
	piecewise_lst *pl;
	Polyhedron *U;
	exvector vars;
	exvector params;
	exvector allvars;
};

static int guarded_qp_bernstein_coefficients(__isl_take isl_set *set,
	__isl_take isl_qpolynomial *qp, void *user)
{
	struct isl_bound *bound = (struct isl_bound *)user;
	unsigned nvar;
	Polyhedron *D;
	struct barvinok_options *options;
	evalue *e;

	options = barvinok_options_new_with_defaults();

	if (!set || !qp)
		goto error;

	nvar = isl_set_dim(set, isl_dim_set);

	e = isl_qpolynomial_to_evalue(qp);
	if (!e)
		goto error;

	set = isl_set_make_disjoint(set);
	D = isl_set_to_polylib(set);

	bound->vars = constructVariableVector(nvar, "v");
	bound->allvars = bound->vars;
	bound->allvars.insert(bound->allvars.end(),
				bound->params.begin(), bound->params.end());

	bound->pl = bernstein_coefficients_relation(bound->pl, D, e, bound->U,
			bound->allvars, bound->vars, bound->params, options);

	Domain_Free(D);
	evalue_free(e);
	isl_set_free(set);
	isl_qpolynomial_free(qp);
	barvinok_options_free(options);

	return 0;
error:
	isl_set_free(set);
	isl_qpolynomial_free(qp);
	barvinok_options_free(options);
	return -1;
}

__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_upper_bound(
	__isl_take isl_pw_qpolynomial *pwqp)
{
	isl_dim *dim = NULL;
	unsigned nvar;
	unsigned nparam;
	struct isl_bound bound;
	struct isl_pw_qpolynomial_fold *pwf;

	if (!pwqp)
		return NULL;

	dim = isl_pw_qpolynomial_get_dim(pwqp);
	nvar = isl_dim_size(dim, isl_dim_set);
	if (nvar == 0) {
		isl_dim_free(dim);
		return isl_pw_qpolynomial_fold_from_pw_qpolynomial(isl_fold_max,
									pwqp);
	}

	nparam = isl_dim_size(dim, isl_dim_param);
	bound.pl = NULL;
	bound.U = Universe_Polyhedron(nparam);
	bound.params = constructVariableVector(nparam, "p");

	if (isl_pw_qpolynomial_foreach_lifted_piece(pwqp,
				    guarded_qp_bernstein_coefficients, &bound))
		goto error;

	bound.pl->maximize();

	dim = isl_dim_drop(dim, isl_dim_set, 0, nvar);

	bound.pl->sign = 1;
	pwf = isl_pw_qpolynomial_fold_from_ginac(dim, bound.pl, bound.params);

	Polyhedron_Free(bound.U);
	delete bound.pl;

	isl_pw_qpolynomial_free(pwqp);

	return pwf;
error:
	isl_pw_qpolynomial_free(pwqp);
	isl_dim_free(dim);
	return NULL;
}
