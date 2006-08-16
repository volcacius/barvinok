#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <polylib/polylibgmp.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    (POL_NO_DUAL | POL_INTEGER)
#else
#define MAXRAYS  600
#endif

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

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

static void scan_poly(Polyhedron *S, int pos, Value *z)
{
    if (!S) {
	int k;
	value_print(stdout, VALUE_FMT, z[1]);
	for (k=2; k <= pos; ++k) {
	    printf(", ");
	    value_print(stdout,VALUE_FMT,z[k]);
	}
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
	    scan_poly(S->next, pos+1, z);
	}
	value_set_si(z[pos+1], 0);
	value_clear(LB);
	value_clear(UB);
	value_clear(tmp);
    }
}

int main(int argc, char **argv)
{
    Polyhedron *A, *U, *S;
    Value *p;
    int i;

    A = Polyhedron_Read();
    U = Universe_Polyhedron(0);
    S = Polyhedron_Scan(A, U, MAXRAYS);

    p = ALLOCN(Value, A->Dimension+2);
    for(i=0;i<=A->Dimension;i++) {
	value_init(p[i]);
	value_set_si(p[i],0);
    }
    value_init(p[i]);
    value_set_si(p[i],1);

    scan_poly(S, 0, p);

    for(i=0;i<=(A->Dimension+1);i++) 
	value_clear(p[i]);
    free(p);
    Domain_Free(S);
    Polyhedron_Free(A);
    Polyhedron_Free(U);
    return 0;
}
