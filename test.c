#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

#ifdef HAVE_SYS_TIMES_H

#include <sys/times.h>

static void time_diff(struct tms *before, struct tms *after)
{
    long ticks = sysconf(_SC_CLK_TCK);
    printf("User: %g; Sys: %g\n", 
	    (0.0 + after->tms_utime - before->tms_utime) / ticks,
	    (0.0 + after->tms_stime - before->tms_stime) / ticks);
}

#else

struct tms {};
static void times(struct tms* time)
{
}
static void time_diff(struct tms *before, struct tms *after)
{
}

#endif

int main()
{
    int i, nbPol, nbVec, func, j;
    Polyhedron *A, *B, *C, *D, *E, *F, *G;
    char s[128];

    nbPol = nbVec = 0;
    fgets(s, 128, stdin);
    while ((*s=='#') ||
	    ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "V %d", &nbVec)<1)) )
	fgets(s, 128, stdin);

    for (i = 0; i < nbPol; ++i) {
	Matrix *M = Matrix_Read();
	A = Constraints2Polyhedron(M, MAXRAYS);
	Matrix_Free(M);
	fgets(s, 128, stdin);
	while ((*s=='#') || (sscanf(s, "F %d", &func)<1))
	    fgets(s, 128, stdin);

	switch(func) {
	case 0: {
	    Value cb, ck;
	    value_init(cb);
	    value_init(ck);
	    fgets(s, 128, stdin);
	    /* workaround for apparent bug in older gmps */
	    *strchr(s, '\n') = '\0';
	    while ((*s=='#') || (value_read(ck, s) != 0)) {
		fgets(s, 128, stdin);
		/* workaround for apparent bug in older gmps */
		*strchr(s, '\n') = '\0';
	    }
	    barvinok_count(A, &cb, MAXRAYS);
	    if (value_ne(cb, ck))
		return -1;
	    value_clear(cb);
	    value_clear(ck);
	    break;
	}
	case 1:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    B = Polyhedron_Polar(A, MAXRAYS);
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    C = Polyhedron_Polar(B, MAXRAYS);
	    Polyhedron_Print(stdout, P_VALUE_FMT, C);
	    Polyhedron_Free(C);
	    Polyhedron_Free(B);
	    break;
	case 2:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    for (j = 0; j < A->NbRays; ++j) {
		B = supporting_cone(A, j);
		Polyhedron_Print(stdout, P_VALUE_FMT, B);
		Polyhedron_Free(B);
	    }
	    break;
	case 3:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    C = B = NULL;
	    barvinok_decompose(A,&B,&C);
	    puts("Pos:");
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    puts("Neg:");
	    Polyhedron_Print(stdout, P_VALUE_FMT, C);
	    Domain_Free(B);
	    Domain_Free(C);
	    break;
	case 4: {
	    Value cm, cb;
	    struct tms tms_before, tms_between, tms_after;
	    value_init(cm);
	    value_init(cb);
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    times(&tms_before);
	    manual_count(A, &cm);
	    times(&tms_between);
	    barvinok_count(A, &cb, 100);
	    times(&tms_after);
	    printf("manual: ");
	    value_print(stdout, P_VALUE_FMT, cm);
	    puts("");
	    time_diff(&tms_before, &tms_between);
	    printf("Barvinok: ");
	    value_print(stdout, P_VALUE_FMT, cb);
	    puts("");
	    time_diff(&tms_between, &tms_after);
	    value_clear(cm);
	    value_clear(cb);
	    break;
	}
	case 5:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    B = triangularize_cone(A, 100);
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    check_triangulization(A, B);
	    Domain_Free(B);
	    break;
	case 6:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    B = remove_equalities(A);
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    Polyhedron_Free(B);
	    break;
	case 7: {
	    Value f;
	    value_init(f);
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    B = Polyhedron_Reduce(A, &f);
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    Polyhedron_Free(B);
	    printf("factor: ");
	    value_print(stdout, P_VALUE_FMT, f);
	    puts("");
	    value_clear(f);
	    break;
	}
	case 8: {
	    Enumeration *en;
	    Matrix *M = Matrix_Read();
	    char **param_name;
	    C = Constraints2Polyhedron(M, MAXRAYS);
	    Matrix_Free(M);
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    Polyhedron_Print(stdout, P_VALUE_FMT, C);
	    en = barvinok_enumerate(A, C, MAXRAYS);
	    param_name = Read_ParamNames(stdin, C->Dimension);
	    print_evalue(stdout, &en->EP, param_name);
	    Enumeration_Free(en);
	    Polyhedron_Free(C);
	}
	case 9:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    Polyhedron_Polarize(A);
	    C = B = NULL;
	    barvinok_decompose(A,&B,&C);
	    for (D = B; D; D = D->next)
		Polyhedron_Polarize(D);
	    for (D = C; D; D = D->next)
		Polyhedron_Polarize(D);
	    puts("Pos:");
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    puts("Neg:");
	    Polyhedron_Print(stdout, P_VALUE_FMT, C);
	    Domain_Free(B);
	    Domain_Free(C);
	    break;
	}
	Polyhedron_Free(A);
    }
    for (i = 0; i < nbVec; ++i) {
	Vector *V = Vector_Read();
	Matrix *M = unimodular_complete(V);
	Matrix_Print(stdout, P_VALUE_FMT, M);
	Matrix_Free(M);
	Vector_Free(V);
    }

    return 0;
}
