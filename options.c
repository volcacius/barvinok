#include <assert.h>
#include <unistd.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "argp.h"
#include "config.h"

#define MAXRAYS    (POL_NO_DUAL | POL_INTEGER)

#define ALLOC(type) (type*)malloc(sizeof(type))

void barvinok_stats_clear(struct barvinok_stats *stats)
{
    memset(stats, 0, sizeof(*stats));
}

void barvinok_stats_print(struct barvinok_stats *stats, FILE *out)
{
    fprintf(out, "Base cones: %ld\n", stats->base_cones);
    if (stats->volume_simplices)
	fprintf(out, "Volume simplices: %ld\n", stats->volume_simplices);
    if (stats->topcom_chambers) {
	fprintf(out, "TOPCOM empty chambers: %ld\n",
		stats->topcom_empty_chambers);
	fprintf(out, "TOPCOM chambers: %ld\n", stats->topcom_chambers);
	fprintf(out, "TOPCOM distinct chambers: %ld\n",
		stats->topcom_distinct_chambers);
    }
    if (stats->gbr_solved_lps)
	fprintf(out, "LPs solved during GBR: %ld\n", stats->gbr_solved_lps);
    if (stats->bernoulli_sums)
	fprintf(out, "Bernoulli sums: %ld\n", stats->bernoulli_sums);
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
    options->count_sample_infinite = 1;
    options->try_Delaunay_triangulation = 0;

#ifdef HAVE_GINAC
    options->bound = BV_BOUND_BERNSTEIN;
#else
    options->bound = BV_BOUND_RANGE;
#endif
    options->chambers = BV_CHAMBERS_POLYLIB;

    options->polynomial_approximation = BV_APPROX_SIGN_NONE;
    options->approximation_method = BV_APPROX_NONE;
    options->scale_flags = 0;
    options->volume_triangulate = BV_VOL_VERTEX;

#ifdef HAVE_LIBGLPK
    options->gbr_lp_solver = BV_GBR_GLPK;
#elif defined HAVE_LIBCDDGMP
    options->gbr_lp_solver = BV_GBR_CDD;
#else
    options->gbr_lp_solver = BV_GBR_PIP;
#endif

#ifdef HAVE_LIBGLPK
    options->lp_solver = BV_LP_GLPK;
#elif defined HAVE_LIBCDDGMP
    options->lp_solver = BV_LP_CDD;
#else
    options->lp_solver = BV_LP_PIP;
#endif

    options->summation = BV_SUM_LAURENT;

    options->bernstein_optimize = BV_BERNSTEIN_NONE;

    options->bernstein_recurse = BV_BERNSTEIN_FACTORS;

    options->integer_hull = BV_HULL_GBR;

    options->verbose = 0;

    options->print_stats = 0;

    options->gbr_only_first = 0;

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
    SCALE_NARROW,
    SCALE_NARROW2,
    SCALE_CHAMBER,
};

const char *scale_opts[] = {
    "fast",
    "slow",
    "narrow",
    "narrow2",
    "chamber",
    NULL
};

static struct argp_option approx_argp_options[] = {
    { "polynomial-approximation", BV_OPT_POLAPPROX, "lower|upper",	1 },
    { "approximation-method", BV_OPT_APPROX,        "scale|drop|volume|bernoulli",	0,
	"method to use in polynomial approximation [default: drop]" },
    { "scale-options",	    BV_OPT_SCALE,
	"fast|slow,narrow|narrow2,chamber",	0 },
    { "volume-triangulation",	    BV_OPT_VOL,	    "lift|vertex|barycenter",    0,
	"type of triangulation to perform in volume computation [default: vertex]" },
    { 0 }
};

static struct argp_option barvinok_argp_options[] = {
    { "index",		    BV_OPT_MAXINDEX,	    "int",		0,
       "maximal index of simple cones in decomposition" },
    { "primal",	    	    BV_OPT_PRIMAL,  	    0,			0 },
    { "table",	    	    BV_OPT_TABLE,  	    0,			0 },
    { "specialization",	    BV_OPT_SPECIALIZATION,  "[bf|df|random|todd]" },
#ifdef HAVE_GINAC
    { "bound",		    BV_OPT_BOUND,	    "bernstein|range",	0,
	"algorithm to use for computing bounds [default: bernstein]" },
#endif
#ifdef POINTS2TRIANGS_PATH
    { "chamber-decomposition", 	BV_OPT_CHAMBERS,    "polylib|topcom",	0,
	"tool to use for chamber decomposition [default: polylib]" },
#endif
    { "gbr",		    BV_OPT_GBR,
#if defined(HAVE_LIBGLPK) && defined(HAVE_LIBCDDGMP)
	"cdd|glpk|pip|pip-dual",
#elif defined(HAVE_LIBGLPK)
	"glpk|pip|pip-dual",
#elif defined(HAVE_LIBCDDGMP)
	"cdd|pip|pip-dual",
#else
	"pip|pip-dual",
#endif
	0,	"lp solver to use for basis reduction "
#ifdef HAVE_LIBGLPK
		"[default: glpk]"
#elif defined HAVE_LIBCDDGMP
		"[default: cdd]"
#else
		"[default: pip]"
#endif
	},
    { "lp",		    BV_OPT_LP,
#if defined(HAVE_LIBGLPK) && defined(HAVE_LIBCDDGMP)
	"cdd|cddf|glpk|pip|polylib",
#elif defined(HAVE_LIBGLPK)
	"glpk|pip|polylib",
#elif defined(HAVE_LIBCDDGMP)
	"cdd|cddf|pip|polylib",
#else
	"pip|polylib",
#endif
	0,	"lp solver to use "
#if defined(HAVE_LIBGLPK)
	"[default: glpk]",
#elif defined(HAVE_LIBCDDGMP)
	"[default: cdd]",
#else
	"[default: pip]",
#endif
	},
    { "summation",	    BV_OPT_SUM,		"box|bernoulli|euler|laurent", 0,
	"[default: laurent]" },
    { "bernstein-recurse",  BV_OPT_RECURSE,    "none|factors|intervals|full",    0,
	"[default: factors]" },
    { "recurse",	    BV_OPT_RECURSE,    	    "",
	OPTION_ALIAS | OPTION_HIDDEN },
    { "integer-hull",	    BV_OPT_HULL,
#ifdef USE_ZSOLVE
	"gbr|hilbert",
#else
	"gbr",
#endif
	0, "[default: gbr]" },
    { "version",	    'V',		    0,			0 },
    { "verbose",    	    'v' },
    { "print-stats",	    BV_OPT_PRINT_STATS,	0,	0 },
    { 0 }
};

static error_t approx_parse_opt(int key, char *arg, struct argp_state *state)
{
    struct barvinok_options *options = state->input;
    char *subopt;

    switch (key) {
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
	if (options->polynomial_approximation == BV_APPROX_SIGN_NONE)
	    options->polynomial_approximation = BV_APPROX_SIGN_APPROX;
	if (!strcmp(arg, "scale"))
	    options->approximation_method = BV_APPROX_SCALE;
	else if (!strcmp(arg, "drop"))
	    options->approximation_method = BV_APPROX_DROP;
	else if (!strcmp(arg, "volume"))
	    options->approximation_method = BV_APPROX_VOLUME;
	else if (!strcmp(arg, "bernoulli"))
	    options->approximation_method = BV_APPROX_BERNOULLI;
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
		options->scale_flags &= ~BV_APPROX_SCALE_NARROW2;
		break;
	    case SCALE_NARROW2:
		options->scale_flags |= BV_APPROX_SCALE_NARROW2;
		options->scale_flags &= ~BV_APPROX_SCALE_NARROW;
		break;
	    case SCALE_CHAMBER:
		options->scale_flags |= BV_APPROX_SCALE_CHAMBER;
		break;
	    default:
		argp_error(state, "unknown suboption '%s'\n", subopt);
	    }
	break;
    case BV_OPT_VOL:
	if (!strcmp(arg, "lift"))
	    options->volume_triangulate = BV_VOL_LIFT;
	else if (!strcmp(arg, "vertex"))
	    options->volume_triangulate = BV_VOL_VERTEX;
	else if (!strcmp(arg, "barycenter"))
	    options->volume_triangulate = BV_VOL_BARYCENTER;
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

static error_t barvinok_parse_opt(int key, char *arg, struct argp_state *state)
{
    struct barvinok_options *options = state->input;
    char *subopt;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options;
	break;
    case 'v':
	options->verbose++;
	break;
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
	else if (!strcmp(arg, "todd"))
	    options->incremental_specialization = BV_SPECIALIZATION_TODD;
	break;
    case BV_OPT_PRIMAL:
	options->primal = 1;
	break;
    case BV_OPT_TABLE:
	options->lookup_table = 1;
	break;
    case BV_OPT_BOUND:
	if (!strcmp(arg, "bernstein"))
	    options->bound = BV_BOUND_BERNSTEIN;
	if (!strcmp(arg, "range"))
	    options->bound = BV_BOUND_RANGE;
	break;
    case BV_OPT_CHAMBERS:
	if (!strcmp(arg, "polylib"))
	    options->chambers = BV_CHAMBERS_POLYLIB;
	if (!strcmp(arg, "topcom"))
	    options->chambers = BV_CHAMBERS_TOPCOM;
	break;
    case BV_OPT_GBR:
	if (!strcmp(arg, "cdd"))
	    options->gbr_lp_solver = BV_GBR_CDD;
	if (!strcmp(arg, "glpk"))
	    options->gbr_lp_solver = BV_GBR_GLPK;
	if (!strcmp(arg, "pip"))
	    options->gbr_lp_solver = BV_GBR_PIP;
	if (!strcmp(arg, "pip-dual"))
	    options->gbr_lp_solver = BV_GBR_PIP_DUAL;
	break;
    case BV_OPT_LP:
	if (!strcmp(arg, "cdd"))
	    options->lp_solver = BV_LP_CDD;
	if (!strcmp(arg, "cddf"))
	    options->lp_solver = BV_LP_CDDF;
	if (!strcmp(arg, "glpk"))
	    options->lp_solver = BV_LP_GLPK;
	if (!strcmp(arg, "pip"))
	    options->lp_solver = BV_LP_PIP;
	if (!strcmp(arg, "polylib"))
	    options->lp_solver = BV_LP_POLYLIB;
	break;
    case BV_OPT_MAXINDEX:
	options->max_index = strtoul(arg, NULL, 0);
	break;
    case BV_OPT_SUM:
	if (!strcmp(arg, "box"))
	    options->summation = BV_SUM_BOX;
	else if (!strcmp(arg, "barvinok"))
	    options->summation = BV_SUM_BOX;
	else if (!strcmp(arg, "euler"))
	    options->summation = BV_SUM_EULER;
	else if (!strcmp(arg, "bernoulli"))
	    options->summation = BV_SUM_BERNOULLI;
	else if (!strcmp(arg, "laurent"))
	    options->summation = BV_SUM_LAURENT;
	else if (!strcmp(arg, "laurent_old"))
	    options->summation = BV_SUM_LAURENT_OLD;
	else
	    argp_error(state, "unknown summation method '%s'\n", arg);
	break;
    case BV_OPT_RECURSE:
	if (!strcmp(arg, "none"))
	    options->bernstein_recurse = 0;
	else if (!strcmp(arg, "factors"))
	    options->bernstein_recurse = BV_BERNSTEIN_FACTORS;
	else if (!strcmp(arg, "intervals"))
	    options->bernstein_recurse = BV_BERNSTEIN_INTERVALS;
	else if (!strcmp(arg, "full"))
	    options->bernstein_recurse =
		    BV_BERNSTEIN_FACTORS | BV_BERNSTEIN_INTERVALS;
	break;
    case BV_OPT_HULL:
	if (!strcmp(arg, "gbr"))
	    options->integer_hull = BV_HULL_GBR;
	else if (!strcmp(arg, "hilbert"))
	    options->integer_hull = BV_HULL_HILBERT;
	break;
    case BV_OPT_PRINT_STATS:
	options->print_stats = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp approx_argp = {
    approx_argp_options, approx_parse_opt, 0, 0
};

static struct argp_child barvinok_children[] = {
    { &approx_argp,    	0,	"polynomial approximation",	BV_GRP_APPROX },
    { 0 }
};

struct argp barvinok_argp = {
    barvinok_argp_options, barvinok_parse_opt, 0, 0, barvinok_children
};
