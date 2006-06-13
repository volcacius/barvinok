#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "config.h"

/* The input of this example program is similar to that of ehrhart_union
 * in the PolyLib distribution, the difference being that the number of
 * polytopes in the union needs to specified explicitly.
 * The input starts with this number, followed by this number of
 * polytopes in combined data and parameter space, a context polytope
 * in parameter space and (optionally) the names of the parameters.
 * All polytopes are in PolyLib notation.
 */

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

#include "config.h"
#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "series",  no_argument,  0,  's' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *C, *D = NULL;
    int i, npol;
    char **param_name;
    char s[128];
    int c, ind = 0;
    int series = 0;

    while ((c = getopt_long(argc, argv, "sV", options, &ind)) != -1) {
	switch (c) {
	case 's':
	    series = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "%d", &npol)<1))
	fgets(s, 128, stdin);

    for (i = 0; i < npol; ++i) {
	Polyhedron *P;
	M = Matrix_Read();
	P = Constraints2Polyhedron(M, MAXRAYS);
	Matrix_Free(M);
	D = DomainConcat(P, D);
    }
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, D);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    if (series) {
	gen_fun *gf;
	gf = barvinok_enumerate_union_series(D, C, MAXRAYS);
	gf->print(C->Dimension, param_name);
	puts("");
	delete gf;
    } else {
	evalue *EP;
	EP = barvinok_enumerate_union(D, C, MAXRAYS);
	print_evalue(stdout, EP, param_name);
	free_evalue_refs(EP);
	free(EP);
    }
    Free_ParamNames(param_name, C->Dimension);
    Domain_Free(D);
    Polyhedron_Free(C);
    return 0;
}
