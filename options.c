#include <unistd.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "argp.h"
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    (POL_NO_DUAL | POL_INTEGER)
#else
#define MAXRAYS  600
#endif

#ifndef HAVE_LIBGLPK
Vector *Polyhedron_Sample(Polyhedron *P, struct barvinok_options *options)
{
    assert(0);
}
#endif

#define ALLOC(type) (type*)malloc(sizeof(type))

void barvinok_stats_clear(struct barvinok_stats *stats)
{
    stats->base_cones = 0;
}

void barvinok_stats_print(struct barvinok_stats *stats, FILE *out)
{
    fprintf(out, "Base cones: %d\n", stats->base_cones);
}

struct barvinok_options *barvinok_options_new_with_defaults()
{
    struct barvinok_options *options = ALLOC(struct barvinok_options);
    if (!options)
	return NULL;

    options->stats = ALLOC(struct barvinok_stats);
    if (!options->stats) {
	free(options);
	return NULL;
    }

    barvinok_stats_clear(options->stats);

    options->LLL_a = 1;
    options->LLL_b = 1;

    options->MaxRays = MAXRAYS;

#ifdef USE_INCREMENTAL_BF
    options->incremental_specialization = 2;
#elif defined USE_INCREMENTAL_DF
    options->incremental_specialization = 1;
#else
    options->incremental_specialization = 0;
#endif
    options->max_index = 1;
    options->primal = 0;
#ifdef USE_MODULO
    options->lookup_table = 0;
#else
    options->lookup_table = 1;
#endif
#ifdef HAVE_LIBGLPK
    options->count_sample_infinite = 1;
#else
    options->count_sample_infinite = 0;
#endif

    options->polynomial_approximation = BV_APPROX_SIGN_NONE;
    options->approximation_method = BV_APPROX_NONE;

#ifdef HAVE_LIBGLPK
    options->gbr_lp_solver = BV_GBR_GLPK;
#elif defined HAVE_LIBCDDGMP
    options->gbr_lp_solver = BV_GBR_CDD;
#else
    options->gbr_lp_solver = BV_GBR_NONE;
#endif

    options->bernstein_optimize = BV_BERNSTEIN_NONE;

    options->bernstein_recurse = BV_BERNSTEIN_FACTORS;

    return options;
}

void barvinok_options_free(struct barvinok_options *options)
{
    free(options->stats);
    free(options);
}

enum {
    SCALE_FAST,
    SCALE_SLOW,
    SCALE_NARROW
};

const char *scale_opts[] = {
    "fast",
    "slow",
    "narrow",
    NULL
};

struct argp_option barvinok_argp_options[] = {
    { "index",		    BV_OPT_MAXINDEX,	    "int",		0,
       "maximal index of simple cones in decomposition" },
    { "primal",	    	    BV_OPT_PRIMAL,  	    0,			0 },
    { "table",	    	    BV_OPT_TABLE,  	    0,			0 },
    { "specialization",	    BV_OPT_SPECIALIZATION,  "[bf|df|random]",	0 },
    { "polynomial-approximation", BV_OPT_POLAPPROX, "lower|upper",	1 },
    { "approximation-method", BV_OPT_APPROX,        "scale|drop",	0,
	"method to use in polynomial approximation [default: drop]" },
    { "scale-options",	    BV_OPT_SCALE,	    "fast|slow,narrow",	0 },
    { "gbr",		    BV_OPT_GBR,    	    "[cdd]",		0,
      "solver to use for basis reduction" },
    { "version",	    'V',		    0,			0 },
    { 0 }
};

error_t barvinok_parse_opt(int key, char *arg, struct argp_state *state)
{
    struct barvinok_options *options = state->input;
    char *subopt;

    switch (key) {
    case 'V':
	printf(barvinok_version());
	exit(0);
    case BV_OPT_SPECIALIZATION:
	if (!strcmp(arg, "bf"))
	    options->incremental_specialization = BV_SPECIALIZATION_BF;
	else if (!strcmp(arg, "df"))
	    options->incremental_specialization = BV_SPECIALIZATION_DF;
	else if (!strcmp(arg, "random"))
	    options->incremental_specialization = BV_SPECIALIZATION_RANDOM;
	break;
    case BV_OPT_PRIMAL:
	options->primal = 1;
	break;
    case BV_OPT_TABLE:
	options->lookup_table = 1;
	break;
    case BV_OPT_GBR:
	if (!strcmp(arg, "cdd"))
	    options->gbr_lp_solver = BV_GBR_CDD;
	break;
    case BV_OPT_MAXINDEX:
	options->max_index = strtoul(arg, NULL, 0);
	break;
    case BV_OPT_POLAPPROX:
	if (!arg) {
	    options->polynomial_approximation = BV_APPROX_SIGN_APPROX;
	    if (options->approximation_method == BV_APPROX_NONE)
		options->approximation_method = BV_APPROX_SCALE;
	} else {
	    if (!strcmp(arg, "lower"))
		options->polynomial_approximation = BV_APPROX_SIGN_LOWER;
	    else if (!strcmp(arg, "upper"))
		options->polynomial_approximation = BV_APPROX_SIGN_UPPER;
	    if (options->approximation_method == BV_APPROX_NONE)
		options->approximation_method = BV_APPROX_DROP;
	}
	break;
    case BV_OPT_APPROX:
	if (!strcmp(arg, "scale"))
	    options->approximation_method = BV_APPROX_SCALE;
	else if (!strcmp(arg, "drop"))
	    options->approximation_method = BV_APPROX_DROP;
	else
	    argp_error(state, "unknown value for --approximation-method option");
	break;
    case BV_OPT_SCALE:
	options->approximation_method = BV_APPROX_SCALE;
	while (*arg != '\0')
	    switch (getsubopt(&arg, scale_opts, &subopt)) {
	    case SCALE_FAST:
		options->scale_flags |= BV_APPROX_SCALE_FAST;
		break;
	    case SCALE_SLOW:
		options->scale_flags &= ~BV_APPROX_SCALE_FAST;
		break;
	    case SCALE_NARROW:
		options->scale_flags |= BV_APPROX_SCALE_NARROW;
		break;
	    default:
		argp_error(state, "unknown suboption '%s'\n", subopt);
	    }
	break;
    case ARGP_KEY_END:
	if (options->polynomial_approximation == BV_APPROX_SIGN_NONE &&
	    options->approximation_method != BV_APPROX_NONE) {
	    fprintf(stderr,
	"no polynomial approximation selected; reseting approximation method\n");
	    options->approximation_method = BV_APPROX_NONE;
	}
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp barvinok_argp = {
    barvinok_argp_options, barvinok_parse_opt, 0, 0
};
