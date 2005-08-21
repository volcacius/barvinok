#include <assert.h>
#include <dlfcn.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>

/* gcc -shared -g collect_polytopes.c -rdynamic -o libcollect.so -ldl -lc -lgmp */

evalue* barvinok_enumerate_e(Polyhedron *P, 
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    static evalue *(*orig)(Polyhedron *, unsigned, unsigned, unsigned) = NULL;
    static char *prefix = NULL;
    evalue *res;
    static int recurse = 0;
    static int counter = 0;

    if (!orig) {
	void *handle = dlopen("libbarvinok.so", RTLD_LAZY);
	assert(handle);
	orig = dlsym(handle, "barvinok_enumerate_e");
	assert(orig);
	dlclose(handle);

	prefix = getenv("POLYTOPE_PREFIX");
	if (prefix)
	    prefix = strdup(prefix);
    }

    if (prefix && !recurse) {
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
	fprintf(f, "\nE %d\nP %d\n", exist, nparam);
	fclose(f);
    }

    recurse = 1;
    res = orig(P, exist, nparam, MaxRays);
    recurse = 0;

    return res;
}
