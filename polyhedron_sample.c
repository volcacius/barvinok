#include <barvinok/options.h>
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
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    M = Matrix_Read();
    P = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);

    sample = Polyhedron_Sample(P, bv_options);
    if (sample) {
	assert(in_domain(P, sample->p));
	Vector_Print(stdout, P_VALUE_FMT, sample);
	Vector_Free(sample);
    }

    Polyhedron_Free(P);
    free(bv_options);
}
