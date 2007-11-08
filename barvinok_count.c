#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "argp.h"
#include "progname.h"

int main(int argc, char **argv)
{
    Value cb;
    Polyhedron *A;
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    set_program_name(argv[0]);
    argp_parse(&barvinok_argp, argc, argv, 0, 0, options);

    A = Polyhedron_Read(options->MaxRays);
    value_init(cb);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    barvinok_count_with_options(A, &cb, options);
    value_print(stdout, P_VALUE_FMT, cb);
    puts("");
    if (options->print_stats)
	barvinok_stats_print(options->stats, stdout);
    value_clear(cb);
    Polyhedron_Free(A);
    barvinok_options_free(options);
    return 0;
}
