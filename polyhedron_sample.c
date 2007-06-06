#include <barvinok/options.h>
#include <barvinok/sample.h>
#include "config.h"
#include "argp.h"
#include "progname.h"

int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *P;
    Vector *sample;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    set_program_name(argv[0]);
    argp_parse(&barvinok_argp, argc, argv, 0, 0, bv_options);

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
    barvinok_options_free(bv_options);
}
