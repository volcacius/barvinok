#include <polylib/polylibgmp.h>
#include <barvinok/sample.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    (POL_NO_DUAL | POL_INTEGER)
#else
#define MAXRAYS  600
#endif

int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *P;
    Vector *sample;

    M = Matrix_Read();
    P = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);

    sample = Polyhedron_Sample(P, MAXRAYS);
    if (sample) {
	Vector_Print(stdout, P_VALUE_FMT, sample);
	Vector_Free(sample);
    }

    Polyhedron_Free(P);
}
