#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "argp.h"

#define PRINT_STATS  	    (BV_OPT_LAST+1)

struct argp_option argp_options[] = {
    { "print-stats",	    PRINT_STATS,  0,	0 },
    { 0 }
};

struct arguments {
    struct barvinok_options *options;
    int print_stats;
};

error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = arguments->options;
	break;
    case PRINT_STATS:
	arguments->print_stats = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

int main(int argc, char **argv)
{
    Value cb;
    Polyhedron *A;
    struct arguments arguments;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    0,	0,  0 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    arguments.print_stats = 0;
    arguments.options = options;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    A = Polyhedron_Read(options->MaxRays);
    value_init(cb);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    barvinok_count_with_options(A, &cb, options);
    value_print(stdout, P_VALUE_FMT, cb);
    puts("");
    if (arguments.print_stats)
	barvinok_stats_print(options->stats, stdout);
    value_clear(cb);
    Polyhedron_Free(A);
    barvinok_options_free(options);
    return 0;
}
