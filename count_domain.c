#include <assert.h>
#include <dlfcn.h>
#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>

/* gcc -shared  count_domain.c -rdynamic -o libcount.so -ldl -lc */

evalue* barvinok_enumerate_ev(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    static evalue *(*orig)(Polyhedron *, Polyhedron *c, unsigned) = NULL;
    evalue *res;
    int c;

    if (!orig) {
	void *handle = dlopen("libbarvinok.so", RTLD_LAZY);
	assert(handle);
	orig = dlsym(handle, "barvinok_enumerate_ev");
	assert(orig);
	dlclose(handle);
    }

    res = orig(P, C, MaxRays);

    c = value_notzero_p(res->d) ? 0
	: res->x.p->type == partition ? res->x.p->size/2 : 1;
    fprintf(stderr, "COUNT: %d %d %d %d\n", P->Dimension, 
	    C->Dimension, c, P->Dimension - P->NbEq/*, domain_size(P)*/);

    return res;
}
