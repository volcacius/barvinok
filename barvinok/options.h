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

    /* lexmin options */
		/* Check for integer points in domain
		 * 0: no
		 * 1: yes
		 */
    int		lexmin_emptiness_check;

    struct barvinok_stats   stats;
};

struct barvinok_options *barvinok_options_new_with_defaults();

#if defined(__cplusplus)
}
#endif

#endif
