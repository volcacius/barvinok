#include <unistd.h>
#include <sys/times.h>
#include <getopt.h>
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
#include <util.h>
#include <barvinok.h>

struct option options[] = {
    { "convert",   no_argument,  0,  'c' },
    { 0, 0, 0, 0 }
};


int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    Enumeration *en;
    char **param_name;
    int c, ind = 0;
    int convert = 0;

    while ((c = getopt_long(argc, argv, "c", options, &ind)) != -1) {
	switch (c) {
	case 'c':
	    convert = 1;
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
    if (convert) {
	Enumeration_mod2table(en, C->Dimension);
	Enumeration_Print(stdout, en, param_name);
    }
    Enumeration_Free(en);
    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    return 0;
}
