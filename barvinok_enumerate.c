#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
#include <util.h>
#include <barvinok.h>

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
    { "size",      no_argument,  0,  's' },
    { 0, 0, 0, 0 }
};
#endif

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    Enumeration *en;
    char **param_name;
    int c, ind = 0;
    int convert = 0;
    int size = 0;

    while ((c = getopt_long(argc, argv, "cs", options, &ind)) != -1) {
	switch (c) {
	case 'c':
	    convert = 1;
	    break;
	case 's':
	    size = 1;
	    break;
	}
    }

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, 600);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, 600);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    en = barvinok_enumerate(A, C, 600);
    Enumeration_Print(stdout, en, param_name);
    if (size)
	printf("\nSize: %d\n", Enumeration_size(en));
    if (convert) {
	Enumeration_mod2table(en, C->Dimension);
	Enumeration_Print(stdout, en, param_name);
	if (size)
	    printf("\nSize: %d\n", Enumeration_size(en));
    }
    Enumeration_Free(en);
    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    return 0;
}
