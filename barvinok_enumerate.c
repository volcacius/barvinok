#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
#include <util.h>
#include <barvinok.h>
#include "config.h"

/* The input of this example program is the same as that of testehrhart
 * in the PolyLib distribution, i.e., a polytope in combined
 * data and parameter space, a context polytope in parameter space
 * and (optionally) the names of the parameters.
 * Both polytopes are in PolyLib notation.
 */

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    0
#else
#define MAXRAYS  600
#endif

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "convert",   no_argument,  0,  'c' },
    { "floor",     no_argument,  0,  'f' },
    { "size",      no_argument,  0,  's' },
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

    while ((c = getopt_long(argc, argv, "fcs", options, &ind)) != -1) {
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
	}
    }

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    EP = barvinok_enumerate_ev(A, C, MAXRAYS);
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
    return 0;
}
