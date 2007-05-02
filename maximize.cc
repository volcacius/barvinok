#include <iostream>
#include <bernstein/bernstein.h>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>
#include "argp.h"
#include "progname.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"

using std::cout;
using std::cerr;
using std::endl;
using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define OPT_VARS  	    (BV_OPT_LAST+1)
#define OPT_SPLIT  	    (BV_OPT_LAST+2)
#define OPT_MIN  	    (BV_OPT_LAST+3)

struct argp_option argp_options[] = {
    { "split",		    OPT_SPLIT,	"int" },
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to maximize" },
    { "verbose",	    'v',  	0,	0, },
    { "minimize",	    OPT_MIN,  	0, 0,	"minimize instead of maximize"},
    { 0 }
};

struct options {
    struct convert_options   convert;
    struct verify_options    verify;
    char* var_list;
    int verbose;
    int split;
    int minimize;
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
	options->verbose = 0;
	options->split = 0;
	options->minimize = 0;
	break;
    case 'v':
	options->verbose = 1;
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    case OPT_SPLIT:
	options->split = atoi(arg);
	break;
    case OPT_MIN:
	options->minimize = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static int check_poly_max(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options);

struct check_poly_max_data : public check_poly_data {
    Polyhedron	    	**S;
    evalue		 *EP;
    piecewise_lst	 *pl;

    check_poly_max_data(Value *z, evalue *EP, piecewise_lst *pl) :
	    		EP(EP), pl(pl) {
	this->z = z;
	this->check = check_poly_max;
    }
};

static void optimum(Polyhedron *S, int pos, const check_poly_max_data *data,
		    Value *opt, bool& found,
		    const struct verify_options *options)
{
    if (!S) {
	Value c;
	value_init(c);
	value_set_double(c, compute_evalue(data->EP, data->z+1)+.25);
	if (!found) {
	    value_assign(*opt, c);
	    found = true;
	} else {
	    if (options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX) {
		if (value_gt(c, *opt))
		    value_assign(*opt, c);
	    } else {
		if (value_lt(c, *opt))
		    value_assign(*opt, c);
	    }
	}
	value_clear(c);
    } else {
	Value LB, UB;
	int ok;
	value_init(LB);
	value_init(UB);
	ok = !(lower_upper_bounds(1+pos, S, data->z, &LB, &UB));
	assert(ok);
	for (; value_le(LB, UB); value_increment(LB, LB)) {
	    value_assign(data->z[1+pos], LB);
	    optimum(S->next, pos+1, data, opt, found, options);
	}
	value_set_si(data->z[1+pos], 0);
	value_clear(LB);
	value_clear(UB);
    }
}

static void optimum(const check_poly_max_data *data, Value *opt,
		    const struct verify_options *options)
{
    bool found = false;
    for (int i = 0; i < data->EP->x.p->size/2; ++i)
	if (!emptyQ2(data->S[i]))
	    optimum(data->S[i], 0, data, opt, found, options);
    assert(found);
}

static int check_poly_max(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options)
{
    int k;
    int ok;
    const check_poly_max_data *max_data;
    max_data = static_cast<const check_poly_max_data *>(data);
    char *minmax;
    Value m, n, d;
    value_init(m);
    value_init(n);
    value_init(d);

    if (options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX)
	minmax = "max";
    else
	minmax = "min";

    max_data->pl->evaluate(nparam, z, &n, &d);
    if (options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX)
	mpz_fdiv_q(m, n, d);
    else
	mpz_cdiv_q(m, n, d);

    if (options->print_all) {
	printf("%s(", minmax);
	value_print(stdout, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    printf(", ");
	    value_print(stdout, VALUE_FMT, z[k]);
	}
	printf(") = ");
	value_print(stdout, VALUE_FMT, n);
	if (value_notone_p(d)) {
	    printf("/");
	    value_print(stdout, VALUE_FMT, d);
	}
	printf(" (");
	value_print(stdout, VALUE_FMT, m);
	printf(")");
    }

    optimum(max_data, &n, options);

    if (options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX)
	ok = value_ge(m, n);
    else
	ok = value_le(m, n);

    if (options->print_all) {
	printf(", %s(EP) = ", minmax);
	value_print(stdout, VALUE_FMT, n);
	printf(". ");
    }

    if (!ok) {
	printf("\n"); 
	fflush(stdout);
	fprintf(stderr, "Error !\n");
	fprintf(stderr, "%s(", minmax);
	value_print(stderr, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    fprintf(stderr,", ");
	    value_print(stderr, VALUE_FMT, z[k]);
	}
	fprintf(stderr, ") should be ");
	if (options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX)
	    fprintf(stderr, "greater");
	else
	    fprintf(stderr, "smaller");
	fprintf(stderr, " than or equal to ");
	value_print(stderr, VALUE_FMT, n);
	fprintf(stderr, ", while pl eval gives ");
	value_print(stderr, VALUE_FMT, m);
	fprintf(stderr, ".\n");
	cerr << *max_data->pl << endl;
    } else if (options->print_all)
	printf("OK.\n");

    value_clear(m);
    value_clear(n);
    value_clear(d);

    return ok;
}

static int verify(Polyhedron *D, piecewise_lst *pl, evalue *EP,
		  unsigned nvar, unsigned nparam, Vector *p,
		  struct verify_options *options)
{
    Polyhedron *CS, *S;
    unsigned MaxRays = options->barvinok->MaxRays;
    assert(value_zero_p(EP->d));
    assert(EP->x.p->type == partition);
    int ok = 1;

    CS = check_poly_context_scan(NULL, &D, D->Dimension, options);

    check_poly_init(D, options);

    if (!(CS && emptyQ2(CS))) {
	check_poly_max_data data(p->p, EP, pl);
	data.S = ALLOCN(Polyhedron *, EP->x.p->size/2);
	for (int i = 0; i < EP->x.p->size/2; ++i) {
	    Polyhedron *A = EVALUE_DOMAIN(EP->x.p->arr[2*i]);
	    data.S[i] = Polyhedron_Scan(A, D, MaxRays & POL_NO_DUAL ? 0 : MaxRays);
	}
	ok = check_poly(CS, &data, nparam, 0, p->p+1+nvar, options);
	for (int i = 0; i < EP->x.p->size/2; ++i)
	    Domain_Free(data.S[i]);
	free(data.S);
    }

    if (!options->print_all)
	printf("\n");

    if (CS) {
	Domain_Free(CS);
	Domain_Free(D);
    }

    return ok;
}

static int verify(piecewise_lst *pl, evalue *EP, unsigned nvar, unsigned nparam,
		  struct verify_options *options)
{
    Vector *p;

    p = Vector_Alloc(nvar+nparam+2);
    value_set_si(p->p[nvar+nparam+1], 1);

    for (int i = 0; i < pl->list.size(); ++i) {
	int ok = verify(pl->list[i].first, pl, EP, nvar, nparam, p, options);
	if (!ok && !options->continue_on_error)
	    break;
    }

    Vector_Free(p);

    return 0;
}

static int optimize(evalue *EP, char **all_vars, unsigned nvar, unsigned nparam,
		    struct options *options)
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
	if (!options->verbose)
	    print_solution = 0;
    }

    if (options->minimize)
	options->verify.barvinok->bernstein_optimize = BV_BERNSTEIN_MIN;
    else
	options->verify.barvinok->bernstein_optimize = BV_BERNSTEIN_MAX;
    pl = evalue_bernstein_coefficients(NULL, EP, U, params,
				       options->verify.barvinok);
    assert(pl);
    if (options->minimize)
    	pl->minimize();
    else
    	pl->maximize();
    if (print_solution)
	cout << *pl << endl;
    if (options->verify.verify)
	result = verify(pl, EP, nvar, nparam, &options->verify);
    delete pl;

    Polyhedron_Free(U);

    return result;
}

int main(int argc, char **argv)
{
    evalue *EP;
    char **all_vars = NULL;
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

    evalue_convert(EP, &options.convert, options.verbose, nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else
	result = optimize(EP, all_vars, nvar, nparam, &options);

    free_evalue_refs(EP);
    free(EP);

    if (options.var_list)
	free(options.var_list);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return result;
}
