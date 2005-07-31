#include <unistd.h>
#include <stdlib.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    0
#else
#define MAXRAYS  600
#endif

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
    Matrix *M;
    int c, ind = 0;

    while ((c = getopt_long(argc, argv, "V", options, &ind)) != -1) {
	switch (c) {
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

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
