#include <unistd.h>
#include <sys/times.h>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
}
#include <util.h>
#include <barvinok.h>
#include <barvinok2.h>
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
    { "explicit",  no_argument,  0,  'e' },
    { 0, 0, 0, 0 }
};
#endif

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    evalue *EP;
    char **param_name;
    gen_fun *gf;
    int c, ind = 0;
    int function = 0;

    while ((c = getopt_long(argc, argv, "e", options, &ind)) != -1) {
	switch (c) {
	case 'e':
	    function = 1;
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
    gf = barvinok_series(A, C, MAXRAYS);
    gf->print(C->Dimension, param_name);
    puts("");
    if (function) {
	EP = *gf;
	print_evalue(stdout, EP, param_name);
    }
    delete gf;
    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    return 0;
}
