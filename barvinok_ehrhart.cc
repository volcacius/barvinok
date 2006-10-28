#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"

/* The input of this example program is a polytope in PolyLib notation,
 * i.e., an n by d+2 matrix of the n constraints A x + b >= 0 defining
 * the polytope * sitting in a d-dimensional space.  The first column
 * is 1 for an inequality and 0 for an equality.  b is placed in the
 * final column.
 * Alternatively, if the matrix is preceded by the word "vertices"
 * on a line by itself, it will be interpreted as a list of vertices
 * in PolyLib notation, i.e., an n by (d+2) matrix, where n is
 * the number of vertices/rays and d the dimension.  The first column is
 * 0 for lines and 1 for vertices/rays.  The final column is the denominator
 * or 0 for rays.  Note that for barvinok_ehrhart, the first column
 * should always be 1.
 */

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "convert",   no_argument,  0,  'c' },
    { "floor",     no_argument,  0,  'f' },
    { "series",    no_argument,  0,  's' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

int main(int argc, char **argv)
{
    Polyhedron *A, *C, *U;
    char **param_name;
    int c, ind = 0;
    int convert = 0;
    int floor = 0;
    int series = 0;
    barvinok_options *bv_options = barvinok_options_new_with_defaults();

    while ((c = getopt_long(argc, argv, "sfcV", options, &ind)) != -1) {
	switch (c) {
	case 's':
	    series = 1;
	    break;
	case 'c':
	    convert = 1;
	    break;
	case 'f':
	    floor = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    A = Polyhedron_Read(bv_options->MaxRays);
    param_name = Read_ParamNames(stdin, 1);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    C = Cone_over_Polyhedron(A);
    U = Universe_Polyhedron(1);
    if (series) {
	gen_fun *gf;
	gf = barvinok_series_with_options(C, U, bv_options);
	gf->print(std::cout, U->Dimension, param_name);
	puts("");
	delete gf;
    } else {
	evalue *EP;
	/* A (conceptually) obvious optimization would be to pass in
	 * the parametric vertices, which are just n times the original
	 * vertices, rather than letting barvinok_enumerate_ev (re)compute
	 * them through Polyhedron2Param_SimplifiedDomain.
	 */
	EP = barvinok_enumerate_with_options(C, U, bv_options);
	print_evalue(stdout, EP, param_name);
	if (floor) {
	    fprintf(stderr, "WARNING: floor conversion not supported\n");
	    evalue_frac2floor(EP);
	    print_evalue(stdout, EP, param_name);
	} else if (convert) {
	    evalue_mod2table(EP, C->Dimension);
	    print_evalue(stdout, EP, param_name);
	}
	free_evalue_refs(EP);
	free(EP);
    }
    Free_ParamNames(param_name, 1);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    Polyhedron_Free(U);
    free(bv_options);
    return 0;
}
