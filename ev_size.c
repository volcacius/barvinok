#include <assert.h>
#include <dlfcn.h>
#include <string.h>
#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>

static long long value_bitsize(Value v);
static long long domain_bitsize(Polyhedron *D);
static long long enode_bitsize(enode *p);
static long long evalue_bitsize(evalue *e);
static long long Enumeration_bitsize(Enumeration *en);

long long value_bitsize(Value v)
{
    int last = (v[0]._mp_size > 0 ? v[0]._mp_size : -v[0]._mp_size) - 1;
    if (last == -1)
	return 0;

    return last * sizeof(v[0]._mp_d[0]) * 8 + ffs(v[0]._mp_d[last]) + 1;
}

long long domain_bitsize(Polyhedron *D)
{
    int i, j;
    long long s = ffs(D->Dimension) + ffs(D->NbConstraints);

    for (i = 0; i < D->NbConstraints; ++i)
	for (j = 0; j < D->Dimension+2; ++j)
	    s += value_bitsize(D->Constraint[i][j]);

    return D->next ? s+domain_bitsize(D->next) : s;
}

long long enode_bitsize(enode *p)
{
    //long long s = (sizeof(*p) - sizeof(p->arr[0])) * 8;
    long long s = ffs(p->size);
    int i;

    if (p->type == partition)
	for (i = 0; i < p->size/2; ++i) {
	    s += domain_bitsize(EVALUE_DOMAIN(p->arr[2*i]));
	    s += evalue_bitsize(&p->arr[2*i+1]);
	}
    else
	for (i = 0; i < p->size; ++i) {
	    s += evalue_bitsize(&p->arr[i]);
	}
    return s;
}

long long evalue_bitsize(evalue *e)
{
    //long long s = sizeof(*e) * 8;
    long long s = 0;
    s += value_bitsize(e->d);
    if (value_notzero_p(e->d))
	s += value_bitsize(e->x.n);
    else
	s += enode_bitsize(e->x.p);
    return s;
}

long long Enumeration_bitsize(Enumeration *en)
{
    long long s = 0;

    for (; en; en = en->next) {
	s += domain_bitsize(en->ValidityDomain);
	s += evalue_bitsize(&en->EP);
    }
    return s;
}

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

    fprintf(stderr, "SIZE: %d %lld %lld\n", P->Dimension - P->NbEq,
	domain_bitsize(P) + domain_bitsize(C), evalue_bitsize(res));

    return res;
}

#undef Enumeration

Enumeration *Polyhedron_Enumerate(Polyhedron *P,Polyhedron *C,unsigned MAXRAYS,char **param_name)
{
    static Enumeration *(*orig)(Polyhedron *, Polyhedron *c, unsigned, char **) = NULL;
    Enumeration *res;
    int c;

    if (!orig) {
	void *handle = dlopen("libpolylibgmp.so", RTLD_LAZY);
	assert(handle);
	orig = dlsym(handle, "Polyhedron_Enumerate");
	assert(orig);
	dlclose(handle);
    }

    res = orig(P, C, MAXRAYS, param_name);

    fprintf(stderr, "SIZE: %d %lld %lld\n", P->Dimension - P->NbEq,
	domain_bitsize(P) + domain_bitsize(C), Enumeration_bitsize(res));

    return res;
}
