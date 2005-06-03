#include <unistd.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    0
#else
#define MAXRAYS  600
#endif

int main()
{
    Value cb;
    Polyhedron *A;
    Matrix *M;

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    value_init(cb);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    barvinok_count(A, &cb, 100);
    value_print(stdout, P_VALUE_FMT, cb);
    puts("");
    value_clear(cb);
    Polyhedron_Free(A);
    return 0;
}
