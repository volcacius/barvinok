#ifndef BARVINOK_OPTIONS_H
#define BARVINOK_OPTIONS_H

#include <stdio.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_stats {
    long	base_cones;
    long	volume_simplices;
    long	topcom_empty_chambers;
    long	topcom_chambers;
    long	topcom_distinct_chambers;
    long	gbr_solved_lps;
    long	bernoulli_sums;
};

void barvinok_stats_clear(struct barvinok_stats *stats);
void barvinok_stats_print(struct barvinok_stats *stats, FILE *out);

struct barvinok_options {
    struct barvinok_stats   *stats;

    /* PolyLib options */
    unsigned	MaxRays;

    /* NTL options */
		/* LLL reduction parameter delta=LLL_a/LLL_b */
    long	LLL_a;
    long	LLL_b;

    /* barvinok options */
    #define	BV_SPECIALIZATION_BF		2
    #define	BV_SPECIALIZATION_DF		1
    #define	BV_SPECIALIZATION_RANDOM	0
    #define	BV_SPECIALIZATION_TODD		3
    int		incremental_specialization;

    unsigned long   	    max_index;
    int		    	    primal;
    int		    	    lookup_table;
    int		    	    count_sample_infinite;

    int			    try_Delaunay_triangulation;

    #define	BV_APPROX_SIGN_NONE	0
    #define	BV_APPROX_SIGN_APPROX	1
    #define	BV_APPROX_SIGN_LOWER	2
    #define	BV_APPROX_SIGN_UPPER	3
    int			    polynomial_approximation;
    #define	BV_APPROX_NONE		0
    #define	BV_APPROX_DROP		1
    #define	BV_APPROX_SCALE		2
    #define	BV_APPROX_VOLUME	3
    #define	BV_APPROX_BERNOULLI	4
    int			    approximation_method;
    #define	BV_APPROX_SCALE_FAST	(1 << 0)
    #define	BV_APPROX_SCALE_NARROW	(1 << 1)
    #define	BV_APPROX_SCALE_NARROW2	(1 << 2)
    #define	BV_APPROX_SCALE_CHAMBER	(1 << 3)
    int			    scale_flags;
    #define	BV_VOL_LIFT		0
    #define	BV_VOL_VERTEX		1
    #define	BV_VOL_BARYCENTER	2
    int			    volume_triangulate;

    /* basis reduction options */
    #define	BV_GBR_PIP	0
    #define	BV_GBR_GLPK	1
    #define	BV_GBR_CDD	2
    #define	BV_GBR_PIP_DUAL	3
    int		gbr_lp_solver;

    /* bernstein options */
    #define	BV_BERNSTEIN_NONE   0
    #define	BV_BERNSTEIN_MAX    1
    #define	BV_BERNSTEIN_MIN   -1
    int		bernstein_optimize;

    #define	BV_BERNSTEIN_FACTORS	1
    #define	BV_BERNSTEIN_INTERVALS	2
    int		bernstein_recurse;

    #define	BV_LP_POLYLIB		0
    #define	BV_LP_GLPK		1
    #define	BV_LP_CDD		2
    #define	BV_LP_CDDF		3
    #define	BV_LP_PIP		4
    int		lp_solver;

    #define	BV_SUM_BOX		0
    /* deprecated */
    #define	BV_SUM_BARVINOK		0
    #define	BV_SUM_EULER		1
    #define	BV_SUM_BERNOULLI	2
    #define	BV_SUM_LAURENT		3
    int		summation;

    #define	BV_CHAMBERS_POLYLIB	0
    #define	BV_CHAMBERS_TOPCOM	1
    int		chambers;

    #define	BV_HULL_GBR		0
    #define	BV_HULL_HILBERT		1
    int		integer_hull;

    int		verbose;

    int		print_stats;

    int		gbr_only_first;
};

struct barvinok_options *barvinok_options_new_with_defaults();
void barvinok_options_free(struct barvinok_options *options);

#define BV_OPT_SPECIALIZATION	256
#define BV_OPT_PRIMAL		257
#define BV_OPT_TABLE		258
#define BV_OPT_GBR		259
#define BV_OPT_MAXINDEX		260
#define BV_OPT_POLAPPROX	261
#define BV_OPT_APPROX		262
#define BV_OPT_SCALE		263
#define BV_OPT_VOL		264
#define BV_OPT_RECURSE		265
#define BV_OPT_LP		266
#define BV_OPT_SUM		267
#define BV_OPT_CHAMBERS		268
#define BV_OPT_HULL		269
#define BV_OPT_PRINT_STATS	270
#define BV_OPT_LAST		270

#define BV_GRP_APPROX		1
#define BV_GRP_LAST		1

struct argp;
extern struct argp barvinok_argp;

#if defined(__cplusplus)
}
#endif

#endif
