#include <assert.h>
#include <dlfcn.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>

/* gcc -shared -g collect_polytopes2.c -rdynamic -o libcollect2.so -ldl -lc -lgmp */

evalue* barvinok_enumerate_ev(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    static evalue *(*orig)(Polyhedron *, Polyhedron *c, unsigned) = NULL;
    static char *prefix = NULL;
    evalue *res;
    static int counter = 0;

    fprintf(stderr, "IN INSTRUMENTED BARVINOK_ENUMERATE_EV");
    fflush(stderr);
    if (!orig) {
	void *handle = dlopen("libbarvinok.so", RTLD_LAZY);
	assert(handle);
	orig = dlsym(handle, "barvinok_enumerate_ev");
	assert(orig);
	dlclose(handle);

	prefix = getenv("POLYTOPE_PREFIX");
	if (prefix)
	    prefix = strdup(prefix);
    }

    if (prefix) {
	char path[PATH_MAX];
	FILE *f;
	int i, j;
	unsigned nr, nc;

	snprintf(path, PATH_MAX, "%s%05d", prefix, counter++);
	f = fopen(path, "w");
	fprintf(f, "%d %d\n", nr=P->NbConstraints, nc=P->Dimension+2);
	for (i=0; i < nr; i++) {
	    for (j=0; j < nc; j++) {
		value_print(f," "P_VALUE_FMT" ", P->Constraint[i][j]);
	    }
	    fprintf(f, "\n");
	}
	fprintf(f, "\n0 %d\n", C->Dimension+2);
	fclose(f);
    }

    res = orig(P, C, MaxRays);

    return res;
}
