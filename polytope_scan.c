#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <barvinok/util.h>
#include <barvinok/options.h>
#include <barvinok/basis_reduction.h>
#include "argp.h"

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

struct argp_option argp_options[] = {
#if defined(HAVE_LIBCDDGMP) || defined(HAVE_LIBGLPK)
    { "direct",		'd',  0,	0,
	"don't apply generalized basis reduction first" },
#endif
    { 0 }
};

struct options {
    struct barvinok_options *barvinok;
    int direct;
};

error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->barvinok;
#if defined(HAVE_LIBCDDGMP) || defined(HAVE_LIBGLPK)
	options->direct = 0;
#else
	options->direct = 1;
#endif
	break;
    case 'd':
	options->direct = 1;
	break;
    case 'V':
	printf(barvinok_version());
	exit(0);
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
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
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    0,	0,  0 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct options options;

    options.barvinok = bv_options;

    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    A = Polyhedron_Read(bv_options->MaxRays);

    if (options.direct) {
	inv = Identity(A->Dimension+1);
	P = A;
    } else {
	basis = Polyhedron_Reduced_Basis(A, bv_options);

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

	P = Polyhedron_Preimage(A, inv, bv_options->MaxRays);
	Polyhedron_Free(A);
    }

    U = Universe_Polyhedron(0);
    S = Polyhedron_Scan(P, U, bv_options->MaxRays);

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
    barvinok_options_free(bv_options);
    return 0;
}
