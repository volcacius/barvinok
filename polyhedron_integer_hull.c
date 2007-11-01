#include <barvinok/polylib.h>
#include <barvinok/options.h>
#include "argp.h"
#include "progname.h"
#include "hull.h"

int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *C;
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    set_program_name(argv[0]);
    argp_parse(&barvinok_argp, argc, argv, 0, 0, options);

    M = Matrix_Read();
    C = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);

    M = Polyhedron_Integer_Hull(C, options);

    Polyhedron_Free(C);

    Matrix_Print(stdout, P_VALUE_FMT, M);
    Matrix_Free(M);

    barvinok_options_free(options);
    return 0;
}
