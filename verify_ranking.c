#include <polylib/polylibgmp.h>
#include <barvinok/util.h>

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

#include "config.h"
#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "continue",  no_argument,  0,  'k' },
    { "max",  no_argument,  0,  'M' },
    { "live",  no_argument,  0,  'l' },
    { "verbose",  no_argument,  0,  'v' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

static int live = 0;
static int print_max = 0;
static int verbose = 0;
static int keep_going = 0;
Value max;

static int check_ranking(Polyhedron *SP, Polyhedron *SD, Enumeration *en,
			 int pos, int nvar, Value *z, Value *count)
{
    int i;
    int ok;
    Value PLB, PUB, DLB, DUB, LB, UB, tmp, c;

    assert(pos < nvar);

    value_init(PLB); value_init(PUB);
    value_init(DLB); value_init(DUB);
    value_init(LB); value_init(UB);
    value_init(tmp);

    if (SP) {
	ok = !(lower_upper_bounds(1+pos, SP, z, &PLB, &PUB));
	assert(ok);
    }
    if (SD) {
	ok = !(lower_upper_bounds(1+pos, SD, z, &DLB, &DUB));
	assert(ok);
    }
    if (!SD || (SP && value_lt(PLB, DLB)))
	value_assign(LB, PLB);
    else
	value_assign(LB, DLB);
    if (!SD || (SP && value_gt(PUB, DUB)))
	value_assign(UB, PUB);
    else
	value_assign(UB, DUB);

    if (pos == nvar-1)
	value_init(c);

    for(value_assign(tmp,LB); value_le(tmp,UB); value_increment(tmp,tmp)) {
	int inP = value_ge(tmp, PLB) && value_le(tmp, PUB);
	int inD = value_ge(tmp, DLB) && value_le(tmp, DUB);
	if (!inP && !inD)
	    continue;

	value_assign(z[pos+1], tmp);
	if (pos < nvar-1)
	    ok &= check_ranking(inP ? SP->next : NULL, 
			        inD ? SD->next : NULL, 
			        en, pos+1, nvar, z, count);

	if (pos == nvar-1 && inD) {
	    Value *ctmp;

	    value_assign(c,*(ctmp=compute_poly(en, z+1)));
	    value_clear(*ctmp);
	    free(ctmp);

	    if (verbose >= 2) {
		printf("EP( ");
		value_print(stdout, VALUE_FMT, z[1]);
		for (i = 2; i <= nvar; ++i) {
		    printf(", ");
		    value_print(stdout, VALUE_FMT, z[i]);
		}
		printf(" ) = ");
		value_print(stdout, VALUE_FMT, c);
		printf("; count = ");
		value_print(stdout, VALUE_FMT, *count);
		printf("\n");
	    }

	    if (value_ne(c, *count)) {
		printf("\n"); 
		fflush(stdout);
		fprintf(stderr,"Error !\n");
		printf("EP( ");
		value_print(stderr, VALUE_FMT, z[1]);
		for (i = 2; i <= nvar; ++i) {
		    printf(", ");
		    value_print(stderr, VALUE_FMT, z[i]);
		}
		printf(" ) = ");
		value_print(stderr, VALUE_FMT, c);
		printf(" but count = ");
		value_print(stderr, VALUE_FMT, *count);
		printf("\n");
		ok = 0;
	    }

	    if (live)
		value_decrement(*count, *count);
	}

	if (!ok && !keep_going)
	    goto end;

	if (pos == nvar-1 && inP) {
	    value_increment(*count, *count);
	    if (value_gt(*count, max))
		value_assign(max, *count);
	}
    }
    value_set_si(z[pos+1], 0);

end:
    if (pos == nvar-1)
	value_clear(c);

    value_clear(PLB); value_clear(PUB);
    value_clear(DLB); value_clear(DUB);
    value_clear(LB); value_clear(UB);
    value_clear(tmp);

    return ok;
}

int main(int argc,char *argv[])
{
    Matrix *M;
    Polyhedron *P, *D, *C;
    Polyhedron *SP, *SD;
    int nb_parms;
    char **param_name = NULL;
    evalue *EP;
    Enumeration *en;
    Vector *z;
    Value count;
    int c, ind = 0;

    while ((c = getopt_long(argc, argv, "klMvV", options, &ind)) != -1) {
	switch (c) {
	case 'k':
	    keep_going = 1;
	    break;
	case 'M':
	    print_max = 1;
	case 'l':
	    live = 1;
	    break;
	case 'v':
	    ++verbose;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    M = Matrix_Read();
    P = Constraints2Polyhedron(M, MAXRAYS);
    assert(P != NULL);
    Matrix_Free(M);
    M = Matrix_Read();
    D = Constraints2Polyhedron(M, MAXRAYS);
    assert(D != NULL);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, MAXRAYS);
    assert(C != NULL);
    Matrix_Free(M);

    nb_parms = D->Dimension;
    param_name = Read_ParamNames(stdin, nb_parms);

    EP = barvinok_ranking_ev(P, D, D->Dimension-C->Dimension, C, MAXRAYS);
    if (live) {
	evalue mone;
	evalue *EC = barvinok_ranking_ev(D, D, D->Dimension-C->Dimension, 
					 C, MAXRAYS);
	if (verbose >= 2) {
	    puts("EP");
	    print_evalue(stdout, EP, param_name);
	    puts("EC");
	    print_evalue(stdout, EC, param_name);
	}
	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
	emul(&mone, EC);
	eadd(EC, EP);
	free_evalue_refs(&mone);
	free_evalue_refs(EC);
	free(EC);
	reduce_evalue(EP);
    }
    if (verbose >= 1) {
	puts("Enumeration");
	print_evalue(stdout, EP, param_name);
    }
    en = partition2enumeration(EP);

    assert(C->Dimension == 0); /* for now */

    /* S = scanning list of polyhedra */
    SP = Polyhedron_Scan(P, C, MAXRAYS);
    SD = Polyhedron_Scan(D, C, MAXRAYS);

    z = Vector_Alloc(P->Dimension+2);
    value_set_si(z->p[P->Dimension+1], 1);

    if (print_max)
	value_init(max);
    value_init(count);
    check_ranking(SP, SD, en, 0, P->Dimension-C->Dimension, z->p, &count);
    value_clear(count);

    if (print_max) {
	printf("max = ");
	value_print(stdout, VALUE_FMT, max);
	printf("\n");
	value_clear(max);
    }

    Enumeration_Free(en);
    Free_ParamNames(param_name, nb_parms);
    Polyhedron_Free(P);
    Polyhedron_Free(D);
    Polyhedron_Free(C);
    Domain_Free(SP);
    Domain_Free(SD);
    Vector_Free(z);
}
