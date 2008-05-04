#include <assert.h>
#include <iostream>
#include <barvinok/barvinok.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"

using std::cout;
using std::cerr;
using std::endl;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define OPT_VARS  	    (BV_OPT_LAST+1)

struct argp_option argp_options[] = {
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to sum" },
    { 0 }
};

struct options {
    struct convert_options   convert;
    struct verify_options    verify;
    char* var_list;
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
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static int check_poly_sum(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options);

struct check_poly_sum_data : public check_EP_data {
    evalue		 *sum;

    check_poly_sum_data(const evalue *EP, evalue *sum) : sum(sum) {
	this->EP = EP;
	this->cp.check = check_poly_sum;
    }
};

static void sum(Polyhedron *S, int pos, const check_poly_sum_data *data,
		evalue *s, const struct verify_options *options)
{
    if (!S) {
	evalue *e = evalue_eval(data->EP, data->cp.z+1);
	eadd(e, s);
	evalue_free(e);
    } else {
	Value LB, UB;
	int ok;
	value_init(LB);
	value_init(UB);
	ok = !(lower_upper_bounds(1+pos, S, data->cp.z, &LB, &UB));
	assert(ok);
	for (; value_le(LB, UB); value_increment(LB, LB)) {
	    value_assign(data->cp.z[1+pos], LB);
	    sum(S->next, pos+1, data, s, options);
	}
	value_set_si(data->cp.z[1+pos], 0);
	value_clear(LB);
	value_clear(UB);
    }
}

static evalue *sum(const check_poly_sum_data *data,
			  const struct verify_options *options)
{
    evalue *s = evalue_zero();
    for (int i = 0; i < data->n_S; ++i)
	if (!emptyQ2(data->S[i]))
	    sum(data->S[i], 0, data, s, options);
    return s;
}

static int check_poly_sum(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options)
{
    const check_poly_sum_data *sum_data;
    sum_data = (const check_poly_sum_data *)data;
    evalue *e, *s;
    int k;
    int ok;

    e = evalue_eval(sum_data->sum, z);
    s = sum(sum_data, options);

    ok = eequal(e, s);

    check_poly_print(ok, nparam, z, s->x.n, s->d, e->x.n, e->d,
		     "sum", "sum(EP)", "summation", options);

    evalue_free(s);
    evalue_free(e);

    return ok;
}

static int verify(evalue *EP, evalue *sum, unsigned nvar, unsigned nparam,
		  struct verify_options *options)
{
    check_poly_sum_data data(EP, sum);
    return !check_EP(&data, nvar, nparam, options);
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

    if (options.verify.verify)
	verify_options_set_range(&options.verify, nvar+nparam);

    evalue_convert(EP, &options.convert, bv_options->verbose, nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else {
	evalue *sum = barvinok_summate(EP, nvar, bv_options);
	assert(sum);
	if (options.verify.verify)
	    result = verify(EP, sum, nvar, nparam, &options.verify);
	else
	    print_evalue(stdout, sum, all_vars+nvar);
	evalue_free(sum);
    }

    evalue_free(EP);

    if (options.var_list)
	free(options.var_list);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return result;
}
