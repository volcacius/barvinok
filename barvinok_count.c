#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
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

static Polyhedron *Polyhedron_Read()
{
    int vertices = 0; 
    unsigned NbRows, NbColumns;
    Matrix *M;
    Polyhedron *P;
    char s[128];

    while (fgets(s, sizeof(s), stdin)) {
	if (*s == '#')
	    continue;
	if (strncasecmp(s, "vertices", sizeof("vertices")-1) == 0)
	    vertices = 1;
	if (sscanf(s, "%u %u", &NbRows, &NbColumns) == 2)
	    break;
    }
    if (feof(stdin))
	return NULL;
    M = Matrix_Alloc(NbRows,NbColumns);
    Matrix_Read_Input(M);
    if (vertices)
	P = Rays2Polyhedron(M, MAXRAYS);
    else
	P = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    return P;
}

int main(int argc, char **argv)
{
    Value cb;
    Polyhedron *A;
    int c, ind = 0;

    while ((c = getopt_long(argc, argv, "V", options, &ind)) != -1) {
	switch (c) {
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    A = Polyhedron_Read();
    value_init(cb);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    barvinok_count(A, &cb, MAXRAYS);
    value_print(stdout, P_VALUE_FMT, cb);
    puts("");
    value_clear(cb);
    Polyhedron_Free(A);
    return 0;
}
