#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/basis_reduction.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    (POL_NO_DUAL | POL_INTEGER)
#else
#define MAXRAYS  600
#endif

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "direct",    no_argument,  0,  'd' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

static Polyhedron *Polyhedron_Read()
{
    int vertices = 0; 
    unsigned NbRows, NbColumns;
    Matrix *M;
    Polyhedron *P;
    char s[128];

    while (fgets(s, sizeof(s), stdin)) {
	if (*s == '#')
	    continue;
	if (strncasecmp(s, "vertices", sizeof("vertices")-1) == 0)
	    vertices = 1;
	if (sscanf(s, "%u %u", &NbRows, &NbColumns) == 2)
	    break;
    }
    if (feof(stdin))
	return NULL;
    M = Matrix_Alloc(NbRows,NbColumns);
    Matrix_Read_Input(M);
    if (vertices)
	P = Rays2Polyhedron(M, MAXRAYS);
    else
	P = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    return P;
}

static void scan_poly(Polyhedron *S, int pos, Value *z, Matrix *T)
{
    if (!S) {
	int k;
	Vector *v;

	v = Vector_Alloc(T->NbRows);
	Matrix_Vector_Product(T, z+1, v->p);
	value_print(stdout, VALUE_FMT, v->p[0]);
	for (k=1; k < pos; ++k) {
	    printf(", ");
	    value_print(stdout,VALUE_FMT, v->p[k]);
	}
	Vector_Free(v);
	printf("\n");
    } else {
	int ok;
	Value LB, UB, tmp;
	value_init(LB);
	value_init(UB);
	value_init(tmp);
	ok = !(lower_upper_bounds(1+pos, S, z, &LB, &UB));
	assert(ok);
	for (value_assign(tmp,LB); value_le(tmp,UB); value_increment(tmp,tmp)) {
	    value_assign(z[pos+1], tmp);
	    scan_poly(S->next, pos+1, z, T);
	}
	value_set_si(z[pos+1], 0);
	value_clear(LB);
	value_clear(UB);
	value_clear(tmp);
    }
}

int main(int argc, char **argv)
{
    Polyhedron *A, *P, *U, *S;
    Value *p;
    int i, j, ok;
    Matrix *basis, *T, *inv;
    int c, ind = 0;
    int direct = 0;

    while ((c = getopt_long(argc, argv, "dV", options, &ind)) != -1) {
	switch (c) {
	case 'd':
	    direct = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    A = Polyhedron_Read();

    if (direct) {
	inv = Identity(A->Dimension+1);
	P = A;
    } else {
	basis = Polyhedron_Reduced_Basis(A);

	T = Matrix_Alloc(A->Dimension+1, A->Dimension+1);
	inv = Matrix_Alloc(A->Dimension+1, A->Dimension+1);
	for (i = 0; i < A->Dimension; ++i)
	    for (j = 0; j < A->Dimension; ++j)
		value_assign(T->p[i][j], basis->p[i][j]);
	value_set_si(T->p[A->Dimension][A->Dimension], 1);
	Matrix_Free(basis);

	ok = Matrix_Inverse(T, inv);
	assert(ok);
	Matrix_Free(T);

	P = Polyhedron_Preimage(A, inv, MAXRAYS);
	Polyhedron_Free(A);
    }

    U = Universe_Polyhedron(0);
    S = Polyhedron_Scan(P, U, MAXRAYS);

    p = ALLOCN(Value, P->Dimension+2);
    for(i=0;i<=P->Dimension;i++) {
	value_init(p[i]);
	value_set_si(p[i],0);
    }
    value_init(p[i]);
    value_set_si(p[i],1);

    scan_poly(S, 0, p, inv);

    Matrix_Free(inv);
    for(i=0;i<=(P->Dimension+1);i++) 
	value_clear(p[i]);
    free(p);
    Domain_Free(S);
    Polyhedron_Free(P);
    Polyhedron_Free(U);
    return 0;
}
