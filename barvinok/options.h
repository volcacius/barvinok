#ifndef BARVINOK_OPTIONS_H
#define BARVINOK_OPTIONS_H

#include <stdio.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_stats {
    long	unimodular_cones;
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
    int		incremental_specialization;

    /* basis reduction options */
    #define	BV_GBR_NONE	0
    #define	BV_GBR_GLPK	1
    #define	BV_GBR_CDD	2
    int		gbr_lp_solver;

    /* lexmin options */
		/* Check for integer points in domain
		 */
    #define	BV_LEXMIN_EMPTINESS_CHECK_NONE		0
    #define	BV_LEXMIN_EMPTINESS_CHECK_SAMPLE	1
    #define	BV_LEXMIN_EMPTINESS_CHECK_COUNT		2
    int		lexmin_emptiness_check;
    int		lexmin_reduce;
    #define	BV_LEXMIN_POLYSIGN_POLYLIB	0
    #define	BV_LEXMIN_POLYSIGN_CDD		1
    #define	BV_LEXMIN_POLYSIGN_CDDF		2
    int		lexmin_polysign;

    struct barvinok_stats   stats;
};

struct barvinok_options *barvinok_options_new_with_defaults();

#if defined(__cplusplus)
}
#endif

#endif
