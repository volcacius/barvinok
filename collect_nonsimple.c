#include <assert.h>
#include <dlfcn.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>

/* gcc -shared -g collect_nonsimple.c -rdynamic -o libcollect.so -ldl -lc -lgmp */

Polyhedron* triangularize_cone(Polyhedron *P, unsigned MaxRays)
{
    static Polyhedron *(*orig)(Polyhedron *, unsigned) = NULL;
    static char *prefix = NULL;
    Polyhedron *res;
    static int recurse = 0;
    static int counter = 0;

    if (!orig) {
	void *handle = dlopen("libbarvinok.so", RTLD_LAZY);
	assert(handle);
	orig = dlsym(handle, "triangularize_cone");
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

	snprintf(path, PATH_MAX, "%s-%05d-%05d", prefix, getpid(), counter++);
	f = fopen(path, "w");
	fprintf(f, "%d %d\n", nr=P->NbRays, nc=P->Dimension+2);
	for (i=0; i < nr; i++) {
	    for (j=0; j < nc; j++) {
		value_print(f," "P_VALUE_FMT" ", P->Ray[i][j]);
	    }
	    fprintf(f, "\n");
	}
	fclose(f);
    }

    recurse = 1;
    res = orig(P, MaxRays);
    recurse = 0;

    return res;
}
