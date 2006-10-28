#include <unistd.h>
#include <stdlib.h>
#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"

/* The input of this example program is the same as that of testehrhart
 * in the PolyLib distribution, i.e., a polytope in combined
 * data and parameter space, a context polytope in parameter space
 * and (optionally) the names of the parameters.
 * Both polytopes are in PolyLib notation.
 */

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "convert",   no_argument,  0,  'c' },
    { "floor",     no_argument,  0,  'f' },
    { "size",      no_argument,  0,  's' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    evalue *EP;
    char **param_name;
    int c, ind = 0;
    int convert = 0;
    int floor = 0;
    int size = 0;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    while ((c = getopt_long(argc, argv, "fcsV", options, &ind)) != -1) {
	switch (c) {
	case 'c':
	    convert = 1;
	    break;
	case 'f':
	    floor = 1;
	    break;
	case 's':
	    size = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    EP = barvinok_enumerate_with_options(A, C, bv_options);
    print_evalue(stdout, EP, param_name);
    if (size)
	printf("\nSize: %d\n", evalue_size(EP));
    if (floor) {
	fprintf(stderr, "WARNING: floor conversion not supported\n");
	evalue_frac2floor(EP);
	print_evalue(stdout, EP, param_name);
    } else if (convert) {
	evalue_mod2table(EP, C->Dimension);
	print_evalue(stdout, EP, param_name);
	if (size)
	    printf("\nSize: %d\n", evalue_size(EP));
    }
    free_evalue_refs(EP);
    free(EP);
    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    free(bv_options);
    return 0;
}
