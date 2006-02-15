#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *C, *D = NULL;
    evalue *EP;
    int i, npol;
    char **param_name;
    char s[128];

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
    EP = barvinok_enumerate_union(D, C, MAXRAYS);
    print_evalue(stdout, EP, param_name);
    free_evalue_refs(EP);
    free(EP);
    Free_ParamNames(param_name, C->Dimension);
    Domain_Free(D);
    Polyhedron_Free(C);
    return 0;
}
