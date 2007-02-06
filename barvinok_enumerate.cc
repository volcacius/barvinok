#include <unistd.h>
#include <stdlib.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "fdstream.h"
#include "argp.h"
#include "verify.h"
#include "verif_ehrhart.h"
#include "remove_equalities.h"

#undef CS   /* for Solaris 10 */

/* The input of this example program is the same as that of testehrhart
 * in the PolyLib distribution, i.e., a polytope in combined
 * data and parameter space, a context polytope in parameter space
 * and (optionally) the names of the parameters.
 * Both polytopes are in PolyLib notation.
 */

#define PRINT_STATS  	    (BV_OPT_LAST+1)

struct argp_option argp_options[] = {
    { "convert",   'c', 0, 0, "convert fractionals to periodics" },
    { "floor",     'f', 0, 0, "convert fractionals to floorings" },
    { "size",      'S' },
    { "series",    's', 0, 0, "compute rational generating function" },
    { "explicit",  'e', 0, 0, "convert rgf to psp" },
    { "verbose",    'v' },
    { "print-stats",	    PRINT_STATS,  0,	0 },
    { 0 }
};

struct arguments {
    struct barvinok_options *barvinok;
    int convert;
    int floor;
    int size;
    int series;
    int function;
    int verbose;
    int print_stats;
    struct verify_options    verify;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *options = (struct arguments*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->barvinok;
	state->child_inputs[1] = &options->verify;
	options->convert = 0;
	options->floor = 0;
	options->size = 0;
	options->series = 0;
	options->function = 0;
	options->verbose = 0;
	options->print_stats = 0;
	break;
    case PRINT_STATS:
	options->print_stats = 1;
	break;
    case 'c':
	options->convert = 1;
	break;
    case 'f':
	options->floor = 1;
	break;
    case 'S':
	options->size = 1;
	break;
    case 'e':
	options->function = 1;
	/* fall through */
    case 's':
	options->series = 1;
	break;
    case 'v':
	options->verbose = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct skewed_gen_fun {
    gen_fun *gf;
    /* maps original space to space in which gf is defined */
    Matrix  *T;
    /* equalities in the original space that need to be satisfied for
     * gf to be valid
     */
    Matrix  *eq;
    /* divisibilities in the original space that need to be satisfied for
     * gf to be valid
     */
    Matrix  *div;

    skewed_gen_fun(gen_fun *gf, Matrix *T, Matrix *eq, Matrix *div) :
		    gf(gf), T(T), eq(eq), div(div) {}
    ~skewed_gen_fun() {
	if (T)
	    Matrix_Free(T);
	if (eq)
	    Matrix_Free(eq);
	if (div)
	    Matrix_Free(div);
	delete gf;
    }

    void print(FILE *out, unsigned int nparam, char **param_name) const;
    operator evalue *() const {
	assert(T == NULL && eq == NULL); /* other cases not supported for now */
	return *gf;
    }
    void coefficient(Value* params, Value* c, barvinok_options *options) const;
};

void skewed_gen_fun::print(FILE *out, unsigned int nparam,
			    char **param_name) const
{
    fdostream os(dup(fileno(out)));
    if (T) {
	fprintf(out, "T:\n");
	Matrix_Print(out, P_VALUE_FMT, T);
    }
    if (eq) {
	fprintf(out, "eq:\n");
	Matrix_Print(out, P_VALUE_FMT, eq);
    }
    if (div) {
	fprintf(out, "div:\n");
	Matrix_Print(out, P_VALUE_FMT, div);
    }
    gf->print(os, nparam, param_name);
}

void skewed_gen_fun::coefficient(Value* params, Value* c,
				 barvinok_options *options) const
{
    if (eq) {
	for (int i = 0; i < eq->NbRows; ++i) {
	    Inner_Product(eq->p[i]+1, params, eq->NbColumns-2, eq->p[i]);
	    if (value_notzero_p(eq->p[i][0])) {
		value_set_si(*c, 0);
		return;
	    }
	}
    }
    if (div) {
	Value tmp;
	value_init(tmp);
	for (int i = 0; i < div->NbRows; ++i) {
	    Inner_Product(div->p[i], params, div->NbColumns-1, &tmp);
	    if (!mpz_divisible_p(tmp, div->p[i][div->NbColumns-1])) {
		value_set_si(*c, 0);
		return;
	    }
	}
	value_clear(tmp);
    }

    ZZ coeff;
    if (!T)
	coeff = gf->coefficient(params, options);
    else {
	Vector *p2 = Vector_Alloc(T->NbRows);
	Matrix_Vector_Product(T, params, p2->p);
	if (value_notone_p(p2->p[T->NbRows-1]))
	    Vector_AntiScale(p2->p, p2->p, p2->p[T->NbRows-1], T->NbRows);
	coeff = gf->coefficient(p2->p, options);
	Vector_Free(p2);
    }

    zz2value(coeff, *c);
}

static int check_series(Polyhedron *S, Polyhedron *CS, skewed_gen_fun *gf,
			int nparam, int pos, Value *z, int print_all,
			barvinok_options *options)
{
    int k;
    Value c, tmp;
    Value LB, UB;

    value_init(c);
    value_init(tmp);
    value_init(LB);
    value_init(UB);

    if (pos == nparam) {
	/* Computes the coefficient */
	gf->coefficient(&z[S->Dimension-nparam+1], &c, options);

	/* if c=0 we may be out of context. */
	/* scanning is useless in this case*/

	if (print_all) {
	    printf("EP( ");
	    value_print(stdout,VALUE_FMT,z[S->Dimension-nparam+1]);
	    for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
		printf(", ");
		value_print(stdout,VALUE_FMT,z[k]);
	    }
	    printf(" ) = ");
	    value_print(stdout,VALUE_FMT,c);
	    printf(" ");
	}

	/* Manually count the number of points */
	count_points(1,S,z,&tmp);
	if (print_all) {
	    printf(", count = ");
	    value_print(stdout, P_VALUE_FMT, tmp);
	    printf(". ");
	}

	if (value_ne(tmp,c)) {
	    printf("\n"); 
	    fflush(stdout);
	    fprintf(stderr,"Error !\n");
	    fprintf(stderr,"EP( ");
	    value_print(stderr,VALUE_FMT,z[S->Dimension-nparam+1]);
	    for (k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
		fprintf(stderr,", ");
		value_print(stderr,VALUE_FMT,z[k]);
	    }
	    fprintf(stderr," ) should be ");
	    value_print(stderr,VALUE_FMT,tmp);
	    fprintf(stderr,", while EP eval gives ");
	    value_print(stderr,VALUE_FMT,c);
	    fprintf(stderr,".\n");
#ifndef DONT_BREAK_ON_ERROR
	    value_clear(c); value_clear(tmp);
	    return 0;
#endif
	} else if (print_all)
	    printf("OK.\n");
    } else {
        int ok = 
	  !(lower_upper_bounds(1+pos, CS, &z[S->Dimension-nparam], &LB, &UB));
        assert(ok);
	for (value_assign(tmp,LB); value_le(tmp,UB); value_increment(tmp,tmp)) {
	    if (!print_all) {
		k = VALUE_TO_INT(tmp);
		if(!pos && !(k%st)) {
		    printf("o");
		    fflush(stdout);
		}
	    }
	    value_assign(z[pos+S->Dimension-nparam+1],tmp);
	    if (!check_series(S, CS->next, gf, nparam, pos+1, z, print_all,
			      options)) {
		value_clear(c); value_clear(tmp);
		value_clear(LB);
		value_clear(UB);
		return(0);
	    }
	}
	value_set_si(z[pos+S->Dimension-nparam+1],0);
    }

    value_clear(c);
    value_clear(tmp);
    value_clear(LB);
    value_clear(UB);
    return 1;
}

static int verify(Polyhedron *P, Polyhedron **C, Enumeration *en, skewed_gen_fun *gf,
		   arguments *options)
{
    Polyhedron *CC, *PP, *CS, *S, *U;
    Matrix *C1, *MM;
    Vector *p;
    int result = 0;

    /******* Compute true context *******/
    CC = align_context(*C, P->Dimension, options->barvinok->MaxRays);
    PP = DomainIntersection(P, CC, options->barvinok->MaxRays);
    Domain_Free(CC);
    C1 = Matrix_Alloc((*C)->Dimension+1, P->Dimension+1);

    for (int i = 0; i < C1->NbRows; i++)
	for (int j = 0; j < C1->NbColumns; j++)
	    if (i == j-P->Dimension+(*C)->Dimension)
		value_set_si(C1->p[i][j], 1);
	    else
		value_set_si(C1->p[i][j], 0);
    CC = Polyhedron_Image(PP, C1, options->barvinok->MaxRays);
    Matrix_Free(C1);
    Domain_Free(PP);
    Domain_Free(*C);
    *C = CC;

    /* Intersect context with range */
    if ((*C)->Dimension > 0) {
	MM = Matrix_Alloc(2*(*C)->Dimension, (*C)->Dimension+2);
	for (int i = 0; i < (*C)->Dimension; ++i) {
	    value_set_si(MM->p[2*i][0], 1);
	    value_set_si(MM->p[2*i][1+i], 1);
	    value_set_si(MM->p[2*i][1+(*C)->Dimension], -options->verify.m);
	    value_set_si(MM->p[2*i+1][0], 1);
	    value_set_si(MM->p[2*i+1][1+i], -1);
	    value_set_si(MM->p[2*i+1][1+(*C)->Dimension], options->verify.M);
	}
	CC = AddConstraints(MM->p[0], 2*(*C)->Dimension, *C,
			    options->barvinok->MaxRays);
	U = Universe_Polyhedron(0);
	CS = Polyhedron_Scan(CC, U, options->barvinok->MaxRays);
	Polyhedron_Free(U);
	Polyhedron_Free(CC);
	Matrix_Free(MM);
    } else
	CS = NULL;

    p = Vector_Alloc(P->Dimension+2);
    value_set_si(p->p[P->Dimension+1], 1);

    /* S = scanning list of polyhedra */
    S = Polyhedron_Scan(P, *C, options->barvinok->MaxRays);

    if (!options->verify.print_all)
	if ((*C)->Dimension > 0) {
	    int d = options->verify.M - options->verify.m;
	    if (d > 80)
		st = 1+d/80;
	    else
		st=1;
	    for (int i = options->verify.m; i <= options->verify.M; i += st)
		printf(".");
	    printf( "\r" );
	    fflush(stdout);
	}

    /******* CHECK NOW *********/
    if (S) {
	if (!options->series || options->function) {
	    if (!check_poly(S, CS, en, (*C)->Dimension, 0, p->p,
			    options->verify.print_all))
		result = -1;
	} else {
	    if (!check_series(S, CS, gf, (*C)->Dimension, 0, p->p,
			      options->verify.print_all, options->barvinok))
		result = -1;
	}
	Domain_Free(S);
    }

    if (result == -1)
	fprintf(stderr,"Check failed !\n");
    
    if (!options->verify.print_all)
	printf( "\n" );
  
    Vector_Free(p);
    if (CS)
	Domain_Free(CS);

    return result;
}

static void unimodular_complete(Matrix *M, int row)
{
    Matrix *H, *Q, *U;
    left_hermite(M, &H, &Q, &U);
    Matrix_Free(H);
    Matrix_Free(U);
    for (int r = row; r < M->NbRows; ++r)
	Vector_Copy(Q->p[r], M->p[r], M->NbColumns);
    Matrix_Free(Q);
}

static skewed_gen_fun *series(Polyhedron *P, Polyhedron* C,
				barvinok_options *options)
{
    Polyhedron *C1, *C2;
    gen_fun *gf;
    Matrix *inv = NULL;
    Matrix *eq = NULL;
    Matrix *div = NULL;
    Polyhedron *PT = P;

    /* Compute true context */
    C1 = Polyhedron_Project(P, C->Dimension);
    C2 = DomainIntersection(C, C1, options->MaxRays);
    Polyhedron_Free(C1);

    POL_ENSURE_VERTICES(C2);
    if (C2->NbBid != 0) {
	Polyhedron *T;
	Matrix *M, *Minv, *M2;
	Matrix *CP;
	if (C2->NbEq || P->NbEq) {
	    /* We remove all equalities to be sure all lines are unit vectors */
	    Polyhedron *CT = C2;
	    remove_all_equalities(&PT, &CT, &CP, NULL, C2->Dimension,
				  options->MaxRays);
	    if (CT != C2) {
		Polyhedron_Free(C2);
		C2 = CT;
	    }
	    if (CP) {
		inv = left_inverse(CP, &eq);
		Matrix_Free(CP);

		int d = 0;
		Value tmp;
		value_init(tmp);
		div = Matrix_Alloc(inv->NbRows-1, inv->NbColumns+1);
		for (int i = 0; i < inv->NbRows-1; ++i) {
		    Vector_Gcd(inv->p[i], inv->NbColumns, &tmp);
		    if (mpz_divisible_p(tmp,
					inv->p[inv->NbRows-1][inv->NbColumns-1]))
			continue;
		    Vector_Copy(inv->p[i], div->p[d], inv->NbColumns);
		    value_assign(div->p[d][inv->NbColumns],
				 inv->p[inv->NbRows-1][inv->NbColumns-1]);
		    ++d;
		}
		value_clear(tmp);

		if (!d) {
		    Matrix_Free(div);
		    div = NULL;
		} else
		    div->NbRows = d;
	    }
	}
	POL_ENSURE_VERTICES(C2);

	/* Since we have "compressed" the parameters (in case there were
	 * any equalities), the result is independent of the coordinates in the
	 * coordinate subspace spanned by the lines.  We can therefore assume
	 * these coordinates are zero and compute the inverse image of the map
	 * from a lower dimensional space that adds zeros in the appropriate
	 * places.
	 */
	M = Matrix_Alloc(C2->Dimension+1, C2->Dimension-C2->NbBid+1);
	int k = 0;
	for (int i = 0; i < C2->NbBid; ++i) {
	    int j = First_Non_Zero(C2->Ray[i]+1, C2->Dimension);
	    assert(First_Non_Zero(C2->Ray[i]+1+j+1, C2->Dimension-j-1) == -1);
	    for ( ; k < j; k++)
		value_set_si(M->p[k+i][k], 1);
	}
	for ( ; k < C2->Dimension-C2->NbBid+1; k++)
	    value_set_si(M->p[k+C2->NbBid][k], 1);
	Minv = Transpose(M);
	M2 = align_matrix(M, PT->Dimension + 1);

	T = PT;
	PT = Polyhedron_Preimage(PT, M2, options->MaxRays);
	if (T != P)
	    Polyhedron_Free(T);

	T = C2;
	C2 = Polyhedron_Preimage(C2, M, options->MaxRays);
	if (T != C)
	    Polyhedron_Free(T);

	Matrix_Free(M);
	Matrix_Free(M2);

	if (inv) {
	    Matrix *T = inv;
	    inv = Matrix_Alloc(Minv->NbRows, T->NbColumns);
	    Matrix_Product(Minv, T, inv);
	    Matrix_Free(T);
	    Matrix_Free(Minv);
	} else
	    inv = Minv;
    }
    gf = barvinok_series_with_options(PT, C2, options);
    Polyhedron_Free(C2);
    if (PT != P)
	Polyhedron_Free(PT);
    return new skewed_gen_fun(gf, inv, eq, div);
}

int main(int argc, char **argv)
{
    Polyhedron *A, *C;
    Matrix *M;
    evalue *EP = NULL;
    Enumeration *en = NULL;
    skewed_gen_fun *gf = NULL;
    char **param_name;
    int print_solution = 1;
    int result = 0;
    struct arguments options;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  		0 },
	{ &verify_argp,    	0,	"verification",	1 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    options.barvinok = bv_options;
    argp_parse(&argp, argc, argv, 0, 0, &options);

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    param_name = Read_ParamNames(stdin, C->Dimension);

    if (options.verify.verify) {
	verify_options_set_range(&options.verify, A);
	if (!options.verbose)
	    print_solution = 0;
    }

    if (print_solution) {
	Polyhedron_Print(stdout, P_VALUE_FMT, A);
	Polyhedron_Print(stdout, P_VALUE_FMT, C);
    }

    if (options.series) {
	gf = series(A, C, bv_options);
	if (print_solution) {
	    gf->print(stdout, C->Dimension, param_name);
	    puts("");
	}
	if (options.function) {
	    EP = *gf;
	    if (print_solution)
		print_evalue(stdout, EP, param_name);
	}
    } else {
	EP = barvinok_enumerate_with_options(A, C, bv_options);
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
	if (options.size)
	    printf("\nSize: %d\n", evalue_size(EP));
	if (options.floor) {
	    fprintf(stderr, "WARNING: floor conversion not supported\n");
	    evalue_frac2floor2(EP, 0);
	    print_evalue(stdout, EP, param_name);
	} else if (options.convert) {
	    evalue_mod2table(EP, C->Dimension);
	    print_evalue(stdout, EP, param_name);
	    if (options.size)
		printf("\nSize: %d\n", evalue_size(EP));
	}
    }

    if (options.verify.verify) {
	if (EP) {
	    en = partition2enumeration(EP);
	    EP = NULL;
	}
	result = verify(A, &C, en, gf, &options);
    }

    if (en)
	Enumeration_Free(en);
    if (gf)
	delete gf;
    if (EP) {
	free_evalue_refs(EP);
	free(EP);
    }

    if (options.print_stats)
	barvinok_stats_print(&options.barvinok->stats, stdout);

    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    barvinok_options_free(bv_options);
    return result;
}
