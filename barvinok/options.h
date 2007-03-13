#ifndef BARVINOK_OPTIONS_H
#define BARVINOK_OPTIONS_H

#include <stdio.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_stats {
    long	base_cones;
};

void barvinok_stats_clear(struct barvinok_stats *stats);
void barvinok_stats_print(struct barvinok_stats *stats, FILE *out);

struct barvinok_options {
    /* PolyLib options */
    unsigned	MaxRays;

    /* NTL options */
		/* LLL reduction parameter delta=LLL_a/LLL_b */
    long	LLL_a;
    long	LLL_b;

    /* barvinok options */
		/*
		 * 0: no
		 * 1: depth first
		 * 2: breadth first
		 */
    #define	BV_SPECIALIZATION_BF		2
    #define	BV_SPECIALIZATION_DF		1
    #define	BV_SPECIALIZATION_RANDOM	0
    int		incremental_specialization;

    unsigned long   	    max_index;
    int		    	    primal;
    int		    	    lookup_table;
    int		    	    count_sample_infinite;

    #define	BV_POLAPPROX_NONE	0
    #define	BV_POLAPPROX_PRE_LOWER	1
    #define	BV_POLAPPROX_PRE_UPPER	2
    #define	BV_POLAPPROX_PRE_APPROX	3
    #define	BV_POLAPPROX_LOWER	4
    #define	BV_POLAPPROX_UPPER	5
    int			    polynomial_approximation;

    /* basis reduction options */
    #define	BV_GBR_NONE	0
    #define	BV_GBR_GLPK	1
    #define	BV_GBR_CDD	2
    int		gbr_lp_solver;

    /* bernstein options */
    #define	BV_BERNSTEIN_NONE   0
    #define	BV_BERNSTEIN_MAX    1
    #define	BV_BERNSTEIN_MIN    2
    int		bernstein_optimize;

    #define	BV_BERNSTEIN_FACTORS	1
    #define	BV_BERNSTEIN_INTERVALS	2
    int		bernstein_recurse;

    struct barvinok_stats   *stats;
};

struct barvinok_options *barvinok_options_new_with_defaults();
void barvinok_options_free(struct barvinok_options *options);

#define BV_OPT_SPECIALIZATION	256
#define BV_OPT_PRIMAL		257
#define BV_OPT_TABLE		258
#define BV_OPT_GBR		259
#define BV_OPT_MAXINDEX		260
#define BV_OPT_POLAPPROX	261
#define BV_OPT_LAST		261

struct argp;
extern struct argp barvinok_argp;

#if defined(__cplusplus)
}
#endif

#endif
