#include <assert.h>
#include <stdio.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>

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
	    decompose(A,&B,&C);
	    puts("Pos:");
	    Polyhedron_Print(stdout, P_VALUE_FMT, B);
	    puts("Neg:");
	    Polyhedron_Print(stdout, P_VALUE_FMT, C);
	    Domain_Free(B);
	    Domain_Free(C);
	    break;
	case 4: {
	    Value c;
	    value_init(c);
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    count(A, &c);
	    value_print(stdout, P_VALUE_FMT, c);
	    puts("");
	    value_clear(c);
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
	}
	Polyhedron_Free(A);
    }

    return 0;
}
