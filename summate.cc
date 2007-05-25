#include <iostream>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"
#include "evalue_convert.h"
#include "evalue_read.h"

using std::cout;
using std::cerr;
using std::endl;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define OPT_VARS  	    (BV_OPT_LAST+1)

struct argp_option argp_options[] = {
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to maximize" },
    { "verbose",	    'v',  	0,	0, },
    { 0 }
};

struct options {
    struct convert_options   convert;
    struct barvinok_options  *barvinok;
    char* var_list;
    int verbose;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = &options->convert;
	state->child_inputs[1] = &options->barvinok;
	options->var_list = NULL;
	options->verbose = 0;
	break;
    case 'v':
	options->verbose = 1;
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
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
	{ &barvinok_argp,    	0,	"barvinok options",	3 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    int result = 0;

    options.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    EP = evalue_read_from_file(stdin, options.var_list, &all_vars,
			       &nvar, &nparam, bv_options->MaxRays);
    assert(EP);

    evalue_convert(EP, &options.convert, options.verbose, nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else {
	evalue *sum = esum(EP, nvar);
	print_evalue(stdout, sum, all_vars+nvar);
	free_evalue_refs(sum);
	free(sum);
    }

    free_evalue_refs(EP);
    free(EP);

    if (options.var_list)
	free(options.var_list);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return result;
}
