#include <polylib/polylibgmp.h>
#include "piputil.h"
#include "ev_operations.h"
#include <barvinok.h>

int main(int argc, char **argv)
{
    Polyhedron *A;
    Matrix *M;
    Polyhedron *D, *P, *N;
    char **param_name;
    int exist, nparam, nvar;
    char s[128];
    evalue sum;

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, 600);
    Matrix_Free(M);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "E %d", &exist)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "P %d", &nparam)<1))
	fgets(s, 128, stdin);

    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    printf("exist: %d, nparam: %d\n", exist, nparam);
    param_name = Read_ParamNames(stdin, nparam);

    nvar = A->Dimension - exist - nparam;
    D = pip_lexmin(A, exist, nparam);

    value_init(sum.d);
    evalue_set_si(&sum, 0, 1);
    for (P = D; P; P = N) {
	N = P->next;
	P->next = 0;
	evalue *EP;
	exist = P->Dimension - nvar - nparam;
	EP = barvinok_enumerate_e(P, exist, nparam, 600);
	print_evalue(stderr, EP, param_name);
	eadd(EP, &sum);
    }
    print_evalue(stderr, &sum, param_name);
}
