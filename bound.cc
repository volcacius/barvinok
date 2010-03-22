#include <assert.h>
#include <iostream>
#include <bernstein/bernstein.h>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include <barvinok/bernstein.h>
#include "argp.h"
#include "progname.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"
#include "range.h"

using std::cout;
using std::cerr;
using std::endl;
using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define OPT_VARS  	    (BV_OPT_LAST+1)
#define OPT_SPLIT  	    (BV_OPT_LAST+2)
#define OPT_LOWER  	    (BV_OPT_LAST+3)
#define OPT_ITERATE  	    (BV_OPT_LAST+4)

struct argp_option argp_options[] = {
    { "split",		    OPT_SPLIT,	"int" },
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to compute a bound" },
    { "lower",	    	    OPT_LOWER, 	0, 0,	"compute lower bound instead of upper bound"},
    { "iterate",	    OPT_ITERATE,  	"int", 1,
	"exact result by iterating over domain (of specified maximal size)"},
    { 0 }
};

struct options {
    struct convert_options   convert;
    struct verify_options    verify;
    char* var_list;
    int split;
    int lower;
    int iterate;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = &options->convert;
	state->child_inputs[1] = &options->verify;
	state->child_inputs[2] = options->verify.barvinok;
	options->var_list = NULL;
	options->split = 0;
	options->lower = 0;
	options->iterate = 0;
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    case OPT_SPLIT:
	options->split = atoi(arg);
	break;
    case OPT_LOWER:
	options->lower = 1;
	break;
    case OPT_ITERATE:
	if (!arg)
	    options->iterate = -1;
	else
	    options->iterate = atoi(arg);
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct verify_point_bound {
	struct verify_point_data vpd;
	isl_pw_qpolynomial *pwqp;
	isl_pw_qpolynomial_fold *pwf;
};

static int verify_point(__isl_take isl_point *pnt, void *user)
{
	int i;
	unsigned nparam;
	struct verify_point_bound *vpb = (struct verify_point_bound *) user;
	isl_int v, n, d, b, t;
	isl_pw_qpolynomial *pwqp;
	isl_qpolynomial *bound;
	isl_qpolynomial *opt;
	const char *minmax;
	int sign;
	int ok;
	int cst;
	FILE *out = vpb->vpd.options->print_all ? stdout : stderr;

	vpb->vpd.n--;

	if (vpb->vpd.options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX) {
		minmax = "max";
		sign = 1;
	} else {
		minmax = "min";
		sign = -1;
	}

	isl_int_init(b);
	isl_int_init(t);
	isl_int_init(v);
	isl_int_init(n);
	isl_int_init(d);

	pwqp = isl_pw_qpolynomial_copy(vpb->pwqp);

	nparam = isl_pw_qpolynomial_dim(pwqp, isl_dim_param);
	for (i = 0; i < nparam; ++i) {
		isl_point_get_coordinate(pnt, isl_dim_param, i, &v);
		pwqp = isl_pw_qpolynomial_fix_dim(pwqp, isl_dim_param, i, v);
	}

	bound = isl_pw_qpolynomial_fold_eval(isl_pw_qpolynomial_fold_copy(vpb->pwf),
						isl_point_copy(pnt));

	if (sign > 0)
		opt = isl_pw_qpolynomial_max(pwqp);
	else
		opt = isl_pw_qpolynomial_min(pwqp);

	cst = isl_qpolynomial_is_cst(opt, &n, &d);
	if (cst != 1)
		goto error;
	if (sign > 0)
		isl_int_fdiv_q(v, n, d);
	else
		isl_int_cdiv_q(v, n, d);

	cst = isl_qpolynomial_is_cst(bound, &n, &d);
	if (cst != 1)
		goto error;
	if (sign > 0)
		isl_int_fdiv_q(b, n, d);
	else
		isl_int_cdiv_q(b, n, d);

	if (sign > 0)
		ok = value_ge(b, v);
	else
		ok = value_le(b, v);

	if (vpb->vpd.options->print_all || !ok) {
		fprintf(out, "%s(", minmax);
		for (i = 0; i < nparam; ++i) {
			if (i)
				fprintf(out, ", ");
			isl_point_get_coordinate(pnt, isl_dim_param, i, &t);
			isl_int_print(out, t, 0);
		}
		fprintf(out, ") = ");
		isl_int_print(out, n, 0);
		if (!isl_int_is_one(d)) {
			fprintf(out, "/");
			isl_int_print(out, d, 0);
			fprintf(out, " (");
			isl_int_print(out, b, 0);
			fprintf(out, ")");
		}
		fprintf(out, ", %s(EP) = ", minmax);
		isl_int_print(out, v, 0);
		if (ok)
			fprintf(out, ". OK\n");
		else
			fprintf(out, ". NOT OK\n");
	} else if ((vpb->vpd.n % vpb->vpd.s) == 0) {
		printf("o");
		fflush(stdout);
	}

	if (0) {
error:
		ok = 0;
	}

	isl_qpolynomial_free(bound);
	isl_qpolynomial_free(opt);
	isl_point_free(pnt);

	isl_int_clear(t);
	isl_int_clear(d);
	isl_int_clear(n);
	isl_int_clear(v);
	isl_int_clear(b);

	if (!ok)
		vpb->vpd.error = 1;

	if (vpb->vpd.options->continue_on_error)
		ok = 1;

	return (vpb->vpd.n >= 1 && ok) ? 0 : -1;
}

static int verify(piecewise_lst *pl, evalue *EP, unsigned nvar,
		  const GiNaC::exvector &params,
		  struct verify_options *options)
{
	struct verify_point_bound vpb = { { options } };
	isl_ctx *ctx = isl_ctx_alloc();
	isl_dim *dim;
	isl_set *context;
	int r;

	dim = isl_dim_set_alloc(ctx, params.size(), 0);
	vpb.pwf = isl_pw_qpolynomial_fold_from_ginac(dim, pl, params);
	dim = isl_dim_set_alloc(ctx, nvar + params.size(), 0);
	vpb.pwqp = isl_pw_qpolynomial_from_evalue(dim, EP);
	vpb.pwqp = isl_pw_qpolynomial_move(vpb.pwqp, isl_dim_set, 0,
						isl_dim_param, 0, nvar);
	context = isl_pw_qpolynomial_fold_domain(
					isl_pw_qpolynomial_fold_copy(vpb.pwf));
	context = verify_context_set_bounds(context, options);

	r = verify_point_data_init(&vpb.vpd, context);

	if (r == 0)
		isl_set_foreach_point(context, verify_point, &vpb);
	if (vpb.vpd.error)
		r = -1;

	isl_set_free(context);
	isl_pw_qpolynomial_free(vpb.pwqp);
	isl_pw_qpolynomial_fold_free(vpb.pwf);

	isl_ctx_free(ctx);

	verify_point_data_fini(&vpb.vpd);

	return r;
}

static piecewise_lst *iterate(evalue *EP, unsigned nvar, Polyhedron *C,
				struct options *options)
{
    assert(C->Dimension == 0);
    Vector *p = Vector_Alloc(nvar+2);
    value_set_si(p->p[nvar+1], 1);
    Value opt;

    value_init(opt);
    struct check_EP_data data;
    data.cp.z = p->p;
    data.cp.check = NULL;

    data.EP = EP;
    check_EP_set_scan(&data, C, options->verify.barvinok->MaxRays);
    int sign = options->verify.barvinok->bernstein_optimize;
    evalue_optimum(&data, &opt, sign);
    check_EP_clear_scan(&data);

    exvector params;
    piecewise_lst *pl = new piecewise_lst(params, sign);
    pl->add_guarded_lst(Polyhedron_Copy(C), lst(value2numeric(opt)));

    value_clear(opt);
    Vector_Free(p);

    return pl;
}

/*
 * Split (partition) EP into a partition with (sub)domains containing
 * size integer points or less and a partition with (sub)domains
 * containing more integer points.
 */
static void split_on_domain_size(evalue *EP, evalue **EP_less, evalue **EP_more,
				 int size, barvinok_options *options)
{
    assert(value_zero_p(EP->d));
    assert(EP->x.p->type == partition);
    assert(EP->x.p->size >= 2);

    struct evalue_section *s_less = new evalue_section[EP->x.p->size/2];
    struct evalue_section *s_more = new evalue_section[EP->x.p->size/2];

    int n_less = 0;
    int n_more = 0;

    Value c;
    value_init(c);

    for (int i = 0; i < EP->x.p->size/2; ++i) {
	Polyhedron *D = EVALUE_DOMAIN(EP->x.p->arr[2*i]);
	Polyhedron *D_less = NULL;
	Polyhedron *D_more = NULL;
	Polyhedron **next_less = &D_less;
	Polyhedron **next_more = &D_more;

	for (Polyhedron *P = D; P; P = P->next) {
	    Polyhedron *next = P->next;
	    P->next = NULL;
	    barvinok_count_with_options(P, &c, options);
	    P->next = next;

	    if (value_zero_p(c))
		continue;

	    if (value_pos_p(c) && value_cmp_si(c, size) <= 0) {
		*next_less = Polyhedron_Copy(P);
		next_less = &(*next_less)->next;
	    } else {
		*next_more = Polyhedron_Copy(P);
		next_more = &(*next_more)->next;
	    }
	}

	if (D_less) {
	    s_less[n_less].D = D_less;
	    s_less[n_less].E = evalue_dup(&EP->x.p->arr[2*i+1]);
	    n_less++;
	} else {
	    s_more[n_more].D = D_more;
	    s_more[n_more].E = evalue_dup(&EP->x.p->arr[2*i+1]);
	    n_more++;
	}
    }

    value_clear(c);

    *EP_less = evalue_from_section_array(s_less, n_less);
    *EP_more = evalue_from_section_array(s_more, n_more);

    delete [] s_less;
    delete [] s_more;
}

static piecewise_lst *optimize(evalue *EP, unsigned nvar, Polyhedron *C,
				const exvector& params,
				struct options *options)
{
    piecewise_lst *pl = NULL;
    if (options->iterate > 0) {
	evalue *EP_less = NULL;
	evalue *EP_more = NULL;
	piecewise_lst *pl_more = NULL;
	split_on_domain_size(EP, &EP_less, &EP_more, options->iterate,
				options->verify.barvinok);
	if (!EVALUE_IS_ZERO(*EP_less)) {
	    options->iterate = -1;
	    pl = optimize(EP_less, nvar, C, params, options);
	}
	if (!EVALUE_IS_ZERO(*EP_more)) {
	    options->iterate = 0;
	    pl_more = optimize(EP_more, nvar, C, params, options);
	}
	evalue_free(EP_less);
	evalue_free(EP_more);
	if (!pl)
	    return pl_more;
	if (!pl_more)
	    return pl;
	pl->combine(*pl_more);
	delete pl_more;
	return pl;
    }
    if (options->iterate)
	pl = iterate(EP, nvar, C, options);
    else if (options->verify.barvinok->bound == BV_BOUND_BERNSTEIN) {
	pl = evalue_bernstein_coefficients(NULL, EP, C, params,
					   options->verify.barvinok);
	if (options->lower)
	    pl->minimize();
	else
	    pl->maximize();
    } else
	pl = evalue_range_propagation(NULL, EP, params,
				      options->verify.barvinok);
    return pl;
}

static int optimize(evalue *EP, const char **all_vars,
		    unsigned nvar, unsigned nparam, struct options *options)
{
    Polyhedron *U;
    piecewise_lst *pl = NULL;
    U = Universe_Polyhedron(nparam);
    int print_solution = 1;
    int result = 0;

    exvector params;
    params = constructParameterVector(all_vars+nvar, nparam);

    if (options->verify.verify) {
	verify_options_set_range(&options->verify, nvar+nparam);
	if (!options->verify.barvinok->verbose)
	    print_solution = 0;
    }

    if (options->lower)
	options->verify.barvinok->bernstein_optimize = BV_BERNSTEIN_MIN;
    else
	options->verify.barvinok->bernstein_optimize = BV_BERNSTEIN_MAX;
    pl = optimize(EP, nvar, U, params, options);
    assert(pl);
    if (print_solution)
	cout << *pl << endl;
    if (options->verify.verify) {
	result = verify(pl, EP, nvar, params, &options->verify);
    }
    delete pl;

    Polyhedron_Free(U);

    return result;
}

int main(int argc, char **argv)
{
    evalue *EP;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct options options;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    static struct argp_child argp_children[] = {
	{ &convert_argp,    	0,	"input conversion",	1 },
	{ &verify_argp,    	0,	"verification",		2 },
	{ &barvinok_argp,    	0,	"barvinok options",	3 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    int result = 0;

    options.verify.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    EP = evalue_read_from_file(stdin, options.var_list, &all_vars,
			       &nvar, &nparam, bv_options->MaxRays);
    assert(EP);

    if (options.split)
	evalue_split_periods(EP, options.split, bv_options->MaxRays);

    evalue_convert(EP, &options.convert, bv_options->verbose,
			nvar+nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else
	result = optimize(EP, all_vars, nvar, nparam, &options);

    evalue_free(EP);

    if (options.var_list)
	free(options.var_list);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return result;
}
