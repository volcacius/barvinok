#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>

int main()
{
    Polyhedron *A, *C;
    Matrix *M;
    Enumeration *en;
    char **param_name;

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, 600);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, 600);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    en = barvinok_enumerate(A, C, 600);
    Enumeration_Print(stdout, en, param_name);
    Enumeration_Free(en);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    return 0;
}
