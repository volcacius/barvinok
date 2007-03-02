#include "conversion.h"
#include "evalue_convert.h"
#include "lattice_point.h"

using std::cout;
using std::cerr;
using std::endl;

static struct argp_option argp_options[] = {
    { "convert",   	    'c', 0, 0, "convert fractionals to periodics" },
    { "combine",   	    'C', 0, 0 },
    { "floor",     	    'f', 0, 0, "convert fractionals to floorings" },
    { "list",   	    'l', 0, 0 },
    { "range-reduction",    'R',    0,	    0 },
    0
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct convert_options *options = (struct convert_options *)state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	options->floor = 0;
	options->convert = 0;
	options->combine = 0;
	options->range = 0;
	options->list = 0;
	break;
    case ARGP_KEY_FINI:
	break;
    case 'f':
	options->floor = 1;
	break;
    case 'c':
	options->convert = 1;
	break;
    case 'C':
	options->combine = 1;
	break;
    case 'l':
	options->list = 1;
	break;
    case 'R':
	options->range = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp convert_argp = {
    argp_options, parse_opt, 0, 0
};

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 : 
	  p->type == flooring ? 1 : 0;
}

static Lattice *extract_lattice(evalue *e, int nparam)
{
    int i;
    Lattice *L;
    Matrix *U;
    Vector *X;
    /* For some mysterious reason, SolveDiophantine expects an extra
     * [0 0 0 1] row in its input matrix.
     */
    Matrix *M = Matrix_Alloc(2, nparam+1+1);
    value_set_si(M->p[1][nparam+1], 1);
    evalue_extract_affine(e, M->p[0], M->p[0]+nparam+1, M->p[0]+nparam);
    /* ignore constant */
    value_set_si(M->p[0][nparam+1], 0);
    SolveDiophantine(M, &U, &X);
    Matrix_Free(M);
    Vector_Free(X);
    L = Matrix_Alloc(nparam+1, nparam+1);
    for (i = 0; i < nparam; ++i)
	Vector_Copy(U->p[i], L->p[i], nparam);
    value_set_si(L->p[nparam][nparam], 1);
    Matrix_Free(U);
    return L;
}

/* Returns a lattice such that the quasi-polynomial e can be represented
 * by a list of polynomials, one for each point in the fundamental
 * parallelepiped of the lattice.
 * If e is a polynomial, then this function returns NULL.
 */
static Lattice *extract_common_lattice(evalue *e, Lattice *L, int nparam)
{
    int i, offset;

    if (value_notzero_p(e->d))
	return L;

    assert(e->x.p->type != partition);

    if (e->x.p->type == fractional) {
	Lattice *L2 = extract_lattice(&e->x.p->arr[0], nparam);
	if (!L)
	    L = L2;
	else {
	    Lattice *L3 = LatticeIntersection(L, L2);
	    Matrix_Free(L);
	    Matrix_Free(L2);
	    L = L3;
	}
    }

    offset = type_offset(e->x.p);
    for (i = e->x.p->size-1; i >= offset; --i)
	L = extract_common_lattice(&e->x.p->arr[i], L, nparam);
    return L;
}

/* Construct an evalue dst from src corresponding to the coset represented
 * by coset, a vector of size number of parameters plus one.
 * The final value in this vector should be 1.
 */
static void evalue_coset(const evalue *src, const Vector *coset, evalue *dst)
{
    if (value_notzero_p(src->d)) {
	value_assign(dst->d, src->d);
	value_init(dst->x.n);
	value_assign(dst->x.n, src->x.n);
	return;
    }

    if (src->x.p->type == fractional) {
	evalue f;
	evalue t;
	Vector *c = Vector_Alloc(coset->Size);
	value_init(f.d);
	value_init(f.x.n);
	evalue_extract_affine(&src->x.p->arr[0], c->p, c->p+c->Size-1, &f.d);
	Inner_Product(coset->p, c->p, c->Size, &f.x.n);
	Vector_Free(c);
	mpz_fdiv_r(f.x.n, f.x.n, f.d);

	evalue_set_si(dst, 0, 1);
	for (int i = src->x.p->size-1; i >= 1; --i) {
	    emul(&f, dst);
	    value_init(t.d);
	    evalue_coset(&src->x.p->arr[i], coset, &t);
	    eadd(&t, dst);
	    free_evalue_refs(&t);
	}
	free_evalue_refs(&f);
	return;
    }

    assert(src->x.p->type == polynomial);
    value_set_si(dst->d, 0);
    dst->x.p = new_enode(src->x.p->type, src->x.p->size, src->x.p->pos);
    for (int i = 0; i < src->x.p->size; ++i)
	evalue_coset(&src->x.p->arr[i], coset, &dst->x.p->arr[i]);
}

static void evalue_print_list_evalue(FILE *out, evalue *e, int nparam,
				     char **params)
{
    Lattice *L;
    L = extract_common_lattice(e, NULL, nparam);
    if (!L)
	print_evalue(out, e, params);
    else {
	Vector *coset = Vector_Alloc(nparam+1);
	value_set_si(coset->p[nparam], 1);
	mat_ZZ RT;
	mat_ZZ R;
	matrix2zz(L, RT, nparam, nparam);
	R = transpose(RT);
	mat_ZZ vertices;
	lattice_point(coset->p, R, vertices, to_ulong(abs(determinant(R))), NULL);
	Matrix_Free(L);
	for (int i = 0; i < vertices.NumRows(); ++i) {
	    evalue t;
	    cout << vertices[i] << endl;
	    zz2values(vertices[i], coset->p);
	    value_init(t.d);
	    evalue_coset(e, coset, &t);
	    print_evalue(stdout, &t, params);
	    free_evalue_refs(&t);
	}
	Vector_Free(coset);
    }
}

static void evalue_print_list(FILE *out, evalue *e, int nparam, char **params)
{
    int i;
    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);

    for (i = 0; i < e->x.p->size/2; i++) {
	Print_Domain(out, EVALUE_DOMAIN(e->x.p->arr[2*i]), params);
	evalue_print_list_evalue(out, &e->x.p->arr[2*i+1], nparam, params);
    }
}

void evalue_convert(evalue *EP, struct convert_options *options, unsigned nparam,
		    char **params)
{
    if (options->combine)
	evalue_combine(EP);
    if (options->range)
	evalue_range_reduction(EP);
    if (params)
	print_evalue(stdout, EP, params);
    if (options->floor) {
	fprintf(stderr, "WARNING: floor conversion not supported\n");
	evalue_frac2floor2(EP, 0);
	if (params)
	    print_evalue(stdout, EP, params);
    } else if (options->list && params) {
	evalue_print_list(stdout, EP, nparam, params);
    } else if (options->convert) {
	evalue_mod2table(EP, nparam);
	if (params)
	    print_evalue(stdout, EP, params);
    }
}
