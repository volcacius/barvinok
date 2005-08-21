#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    0
#else
#define MAXRAYS  600
#endif

void dump_polytope(Polyhedron *P)
{
    int i, j;
    unsigned nr, nc;

    fprintf(stdout, "%d %d\n", nr=P->NbConstraints, nc=P->Dimension+2);
    for (i=0; i < nr; i++) {
	for (j=0; j < nc; j++) {
	    value_print(stdout," "P_VALUE_FMT" ", P->Constraint[i][j]);
	}
	fprintf(stdout, "\n");
    }
}

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    Enumeration *en;
    char **param_name;
    int i;

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    param_name = Read_ParamNames(stdin, C->Dimension);
    A = remove_equalities_p(A, A->Dimension-C->Dimension, 0);
    dump_polytope(A);
    puts("");
    dump_polytope(C);
    puts("");
    for (i = 0; i < C->Dimension; ++i)
	printf("%s ", param_name[i]);
    puts("");
    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    return 0;
}
