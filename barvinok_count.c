#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

int main(int argc, char **argv)
{
    Value cb;
    Polyhedron *A;
    int c, ind = 0;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    while ((c = getopt_long(argc, argv, "V", options, &ind)) != -1) {
	switch (c) {
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    A = Polyhedron_Read(bv_options->MaxRays);
    value_init(cb);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    barvinok_count_with_options(A, &cb, bv_options);
    value_print(stdout, P_VALUE_FMT, cb);
    puts("");
    value_clear(cb);
    Polyhedron_Free(A);
    free(bv_options);
    return 0;
}
