#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>

static void time_diff(struct tms *before, struct tms *after)
{
    long ticks = sysconf(_SC_CLK_TCK);
    printf("User: %g; Sys: %g\n", 
	    (0.0 + after->tms_utime - before->tms_utime) / ticks,
	    (0.0 + after->tms_stime - before->tms_stime) / ticks);
}

int main()
{
    int i, nbPol, nbMat, func, j;
    Polyhedron *A, *B, *C, *D, *E, *F, *G;
    char s[128];

    nbPol = nbMat = 0;
    fgets(s, 128, stdin);
    while ((*s=='#') ||
	    ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
	fgets(s, 128, stdin);

    for (i = 0; i < nbPol; ++i) {
	Matrix *M = Matrix_Read();
	A = Constraints2Polyhedron(M, 600);
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
	    while ((*s=='#') || (value_read(ck, s) != 0))
		fgets(s, 128, stdin);
	    barvinok_count(A, &cb);
	    if (value_ne(cb, ck))
		return -1;
	    value_clear(cb);
	    value_clear(ck);
	    break;
	}
	case 1:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    B = Polyhedron_Polar(A, 600);
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    C = Polyhedron_Polar(B, 600);
	    Polyhedron_Print(stdout, P_VALUE_FMT, C);
	    Polyhedron_Free(C);
	    Polyhedron_Free(B);
	    break;
	case 2:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    for (j = 0; j < A->NbRays; ++j) {
		B = supporting_cone(A, j, 600);
		Polyhedron_Print(stdout, P_VALUE_FMT, B);
		Polyhedron_Free(B);
	    }
	    break;
	case 3:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
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
	    barvinok_count(A, &cb);
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
	    B = triangularize_cone(A, 600);
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    for (C = B; C; C = C->next)
		for (D = C->next; D; D = D->next) {
		    F = C->next;
		    G = D->next;
		    C->next = NULL;
		    D->next = NULL;
		    E = DomainIntersection(C, D, 600);
		    Polyhedron_Print(stdout, P_VALUE_FMT, E);
		    assert(E->NbRays == 0 || E->NbEq == 1);
		    Polyhedron_Free(E);
		    C->next = F;
		    D->next = G;
		}
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
	}
	Polyhedron_Free(A);
    }

    return 0;
}
