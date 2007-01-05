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
    stats->unimodular_cones = 0;
}

void barvinok_stats_print(struct barvinok_stats *stats, FILE *out)
{
    fprintf(out, "Unimodular cones: %d\n", stats->unimodular_cones);
}

struct barvinok_options *barvinok_options_new_with_defaults()
{
    struct barvinok_options *options = ALLOC(struct barvinok_options);
    if (!options)
	return NULL;

    barvinok_stats_clear(&options->stats);

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

#ifdef HAVE_LIBGLPK
    options->gbr_lp_solver = BV_GBR_GLPK;
#elif defined HAVE_LIBCDDGMP
    options->gbr_lp_solver = BV_GBR_CDD;
#else
    options->gbr_lp_solver = BV_GBR_NONE;
#endif

    options->lexmin_emptiness_check = BV_LEXMIN_EMPTINESS_CHECK_SAMPLE;
    options->lexmin_reduce = 1;
    options->lexmin_polysign = BV_LEXMIN_POLYSIGN_POLYLIB;

    return options;
}

struct argp_option barvinok_argp_options[] = {
    { "primal",	    	    BV_OPT_PRIMAL,  	    0,			0 },
    { "table",	    	    BV_OPT_TABLE,  	    0,			0 },
    { "specialization",	    BV_OPT_SPECIALIZATION,  "[bf|df|random]",	0 },
    { "gbr",		    BV_OPT_GBR,    	    "[cdd]",		0,
      "solver to use for basis reduction" },
    { "version",	    'V',		    0,			0 },
    { 0 }
};

error_t barvinok_parse_opt(int key, char *arg, struct argp_state *state)
{
    struct barvinok_options *options = state->input;

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
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp barvinok_argp = {
    barvinok_argp_options, barvinok_parse_opt, 0, 0
};
