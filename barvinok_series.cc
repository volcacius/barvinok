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

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    Enumeration *en;
    char **param_name;
    gen_fun *gf;

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
    gf->print();
    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    return 0;
}
