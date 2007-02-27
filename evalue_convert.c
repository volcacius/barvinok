#include "evalue_convert.h"

static struct argp_option argp_options[] = {
    { "convert",   	    'c', 0, 0, "convert fractionals to periodics" },
    { "combine",   	    'C', 0, 0 },
    { "floor",     	    'f', 0, 0, "convert fractionals to floorings" },
    { "range-reduction",    'R',    0,	    0 },
    0
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct convert_options *options = state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	options->floor = 0;
	options->convert = 0;
	options->combine = 0;
	options->range = 0;
	break;
    case ARGP_KEY_FINI:
	break;
    case 'f':
	options->floor = 1;
	break;
    case 'c':
	options->convert = 1;
	break;
    case 'C':
	options->combine = 1;
	break;
    case 'R':
	options->range = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp convert_argp = {
    argp_options, parse_opt, 0, 0
};

void evalue_convert(evalue *EP, struct convert_options *options, unsigned nparam,
		    char **params)
{
    if (options->combine)
	evalue_combine(EP);
    if (options->range)
	evalue_range_reduction(EP);
    if (params)
	print_evalue(stdout, EP, params);
    if (options->floor) {
	fprintf(stderr, "WARNING: floor conversion not supported\n");
	evalue_frac2floor2(EP, 0);
	if (params)
	    print_evalue(stdout, EP, params);
    } else if (options->convert) {
	evalue_mod2table(EP, nparam);
	if (params)
	    print_evalue(stdout, EP, params);
    }
}
