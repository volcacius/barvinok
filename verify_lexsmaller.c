#include <assert.h>
#include <barvinok/util.h>
#include "config.h"

#define MAXRAYS    POL_NO_DUAL

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
static Value max;

#define LS_OK	    1
#define LS_P	    2	    /* continue searching P */
#define LS_D	    4	    /* continue searching D */

static int check_lexsmaller(Polyhedron *SP, Polyhedron *SD, Enumeration *en,
			 int pos, int nvar, Value *zP, Value *zD, Value *zE,
			 Value *count)
{
    int i;
    int ok;
    Value PLB, PUB, DLB, DUB, LB, UB, tmp, c;

    if (!SP && !SD)
	return LS_OK | LS_P | LS_D;

    value_init(PLB); value_init(PUB);
    value_init(DLB); value_init(DUB);
    value_init(LB); value_init(UB);
    value_init(tmp);

    if (SP) {
	ok = !(lower_upper_bounds(1+pos, SP, zP, &PLB, &PUB));
	assert(ok);
    }
    if (SD) {
	ok = !(lower_upper_bounds(1+pos, SD, zD, &DLB, &DUB));
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

    if (SD && !SD->next)
	value_init(c);

    ok = LS_OK | LS_P | LS_D;

    for(value_assign(tmp,LB); value_le(tmp,UB); value_increment(tmp,tmp)) {
	int inP = SP && value_ge(tmp, PLB) && value_le(tmp, PUB);
	int inD = SD && value_ge(tmp, DLB) && value_le(tmp, DUB);
	if (!inP && !inD)
	    continue;

	if (inP)
	    value_assign(zP[pos+1], tmp);
	if (inD)
	    value_assign(zD[pos+1], tmp);
	if (inD && pos < nvar)
	    value_assign(zE[pos], tmp);

	if (inD && !SD->next) {
	    Value *ctmp;

	    value_assign(c,*(ctmp=compute_poly(en, zE)));
	    value_clear(*ctmp);
	    free(ctmp);

	    if (verbose >= 2) {
		printf("EP( ");
		value_print(stdout, VALUE_FMT, zE[0]);
		for (i = 1; i < nvar; ++i) {
		    printf(", ");
		    value_print(stdout, VALUE_FMT, zE[i]);
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
		fprintf(stderr, "EP( ");
		value_print(stderr, VALUE_FMT, zE[0]);
		for (i = 1; i < nvar; ++i) {
		    fprintf(stderr, ", ");
		    value_print(stderr, VALUE_FMT, zE[i]);
		}
		fprintf(stderr, " ) = ");
		value_print(stderr, VALUE_FMT, c);
		fprintf(stderr, " but count = ");
		value_print(stderr, VALUE_FMT, *count);
		printf("\n");
		ok = 0;
	    }

	    if (live)
		value_decrement(*count, *count);

	    ok &= ~LS_D;
	}

	if (pos < nvar-1)
	    ok &= check_lexsmaller(inP ? SP->next : NULL, 
				   inD ? SD->next : NULL, 
				   en, pos+1, nvar, zP, zD, zE, count);
	else {
	    ok &= check_lexsmaller(NULL, inD ? SD->next : NULL, 
				   en, pos+1, nvar, zP, zD, zE, count)
	       &  check_lexsmaller(inP ? SP->next : NULL, NULL, 
				   en, pos+1, nvar, zP, zD, zE, count);
	    if (pos >= nvar && !(ok & LS_D))
		break;
	    if (pos >= nvar && !(ok & LS_P))
		break;
	}

	if (!ok && !keep_going)
	    goto end;

	if (inP && !SP->next) {
	    value_increment(*count, *count);
	    if (value_gt(*count, max))
		value_assign(max, *count);
	    ok &= ~LS_P;
	}
    }
    if (SP)
	value_set_si(zP[pos+1], 0);
    if (SD)
	value_set_si(zD[pos+1], 0);

end:
    if (SD && !SD->next)
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
    const char **param_name = NULL;
    evalue *EP;
    Enumeration *en;
    Vector *zP, *zD, *zE;
    Value count;
    int c, ind = 0;
    char s[128];
    unsigned dim;

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

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "D %u", &dim)<1))
	fgets(s, 128, stdin);

    M = Matrix_Read();
    C = Constraints2Polyhedron(M, MAXRAYS);
    assert(C != NULL);
    Matrix_Free(M);

    nb_parms = D->Dimension;
    param_name = Read_ParamNames(stdin, nb_parms);

    EP = barvinok_lexsmaller_ev(P, D, dim, C, MAXRAYS);
    if (live) {
	evalue mone;
	evalue *EC = barvinok_lexsmaller_ev(D, D, dim, C, MAXRAYS);
	if (verbose >= 2) {
	    puts("EP");
	    print_evalue(stdout, EP, (const char **)param_name);
	    puts("EC");
	    print_evalue(stdout, EC, (const char **)param_name);
	}
	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
	emul(&mone, EC);
	eadd(EC, EP);
	free_evalue_refs(&mone);
	evalue_free(EC);
	reduce_evalue(EP);
    }
    if (verbose >= 1) {
	puts("Enumeration");
	print_evalue(stdout, EP, (const char **)param_name);
    }
    en = partition2enumeration(EP);

    assert(C->Dimension == 0); /* for now */

    /* S = scanning list of polyhedra */
    SP = Polyhedron_Scan(P, C, MAXRAYS);
    SD = Polyhedron_Scan(D, C, MAXRAYS);

    zP = Vector_Alloc(1+P->Dimension+1);
    value_set_si(zP->p[1+P->Dimension], 1);
    zD = Vector_Alloc(1+D->Dimension+1);
    value_set_si(zD->p[1+D->Dimension], 1);
    zE = Vector_Alloc(dim+C->Dimension);

    if (print_max)
	value_init(max);
    value_init(count);
    check_lexsmaller(SP, SD, en, 0, dim, zP->p, zD->p, zE->p, &count);
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
    Vector_Free(zP);
    Vector_Free(zD);
    Vector_Free(zE);
}
