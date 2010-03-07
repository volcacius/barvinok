#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <isl_set_polylib.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "argp.h"
#include "progname.h"
#include "verify.h"
#include "verif_ehrhart.h"
#include "verify_series.h"
#include "remove_equalities.h"
#include "evalue_convert.h"
#include "conversion.h"
#include "skewed_genfun.h"

#undef CS   /* for Solaris 10 */

using std::cout;
using std::endl;

/* The input of this example program is the same as that of testehrhart
 * in the PolyLib distribution, i.e., a polytope in combined
 * data and parameter space, a context polytope in parameter space
 * and (optionally) the names of the parameters.
 * Both polytopes are in PolyLib notation.
 */

struct argp_option argp_options[] = {
    { "size",      'S' },
    { "series",    's', 0, 0, "compute rational generating function" },
    { "explicit",  'e', 0, 0, "convert rgf to psp" },
    { 0 }
};

struct arguments {
    int size;
    int series;
    int function;
    struct verify_options    verify;
    struct convert_options   convert;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *options = (struct arguments*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->verify.barvinok;
	state->child_inputs[1] = &options->verify;
	state->child_inputs[2] = &options->convert;
	options->size = 0;
	options->series = 0;
	options->function = 0;
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
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static __isl_give isl_set *set_bounds(__isl_take isl_set *set,
	const struct verify_options *options)
{
	int i;
	unsigned nparam;
	isl_point *pt, *pt2;
	isl_set *box;

	nparam = isl_set_dim(set, isl_dim_param);

	if (options->r > 0) {
		pt = isl_set_sample_point(isl_set_copy(set));
		pt2 = isl_point_copy(pt);

		for (i = 0; i < nparam; ++i) {
			pt = isl_point_add_ui(pt, isl_dim_param, i, options->r);
			pt2 = isl_point_sub_ui(pt2, isl_dim_param, i, options->r);
		}
	} else {
		isl_int v;

		isl_int_init(v);
		pt = isl_point_zero(isl_set_get_dim(set));
		isl_int_set_si(v, options->m);
		for (i = 0; i < nparam; ++i) 
			pt = isl_point_set_coordinate(pt, isl_dim_param, i, v);

		pt2 = isl_point_zero(isl_set_get_dim(set));
		isl_int_set_si(v, options->M);
		for (i = 0; i < nparam; ++i) 
			pt2 = isl_point_set_coordinate(pt2, isl_dim_param, i, v);

		isl_int_clear(v);
	}

	box = isl_set_box_from_points(pt, pt2);

	return isl_set_intersect(set, box);
}

struct verify_point_data {
	const struct verify_options *options;
	isl_set *set;
	isl_pw_qpolynomial *pwqp;
	int n;
	int s;
	int error;
};

static int verify_point(__isl_take isl_point *pnt, void *user)
{
	struct verify_point_data *vpd = (struct verify_point_data *) user;
	isl_set *set;
	int i;
	unsigned nparam;
	isl_int v, n, d;
	isl_qpolynomial *cnt = NULL;
	int pa = vpd->options->barvinok->polynomial_approximation;
	int cst;
	int ok;
	FILE *out = vpd->options->print_all ? stdout : stderr;

	vpd->n--;

	isl_int_init(v);
	isl_int_init(n);
	isl_int_init(d);
	set = isl_set_copy(vpd->set);
	nparam = isl_set_dim(set, isl_dim_param);
	for (i = 0; i < nparam; ++i) {
		isl_point_get_coordinate(pnt, isl_dim_param, i, &v);
		set = isl_set_fix(set, isl_dim_param, i, v);
	}

	if (isl_set_count(set, &v) < 0)
		goto error;

	cnt = isl_pw_qpolynomial_eval(isl_pw_qpolynomial_copy(vpd->pwqp),
					isl_point_copy(pnt));

	cst = isl_qpolynomial_is_cst(cnt, &n, &d);
	if (cst != 1)
		goto error;

	if (pa == BV_APPROX_SIGN_LOWER)
		isl_int_cdiv_q(n, n, d);
	else if (pa == BV_APPROX_SIGN_UPPER)
		isl_int_fdiv_q(n, n, d);
	else
		isl_int_tdiv_q(n, n, d);

	if (pa == BV_APPROX_SIGN_APPROX)
		/* just accept everything */
		ok = 1;
	else if (pa == BV_APPROX_SIGN_LOWER)
		ok = isl_int_le(n, v);
	else if (pa == BV_APPROX_SIGN_UPPER)
		ok = isl_int_ge(n, v);
	else
		ok = isl_int_eq(n, v);

	if (vpd->options->print_all || !ok) {
		fprintf(out, "EP(");
		for (i = 0; i < nparam; ++i) {
			if (i)
				fprintf(out, ", ");
			isl_point_get_coordinate(pnt, isl_dim_param, i, &d);
			isl_int_print(out, d, 0);
		}
		fprintf(out, ") = ");
		isl_int_print(out, n, 0);
		fprintf(out, ", count = ");
		isl_int_print(out, v, 0);
		if (ok)
			fprintf(out, ". OK\n");
		else
			fprintf(out, ". NOT OK\n");
	} else if ((vpd->n % vpd->s) == 0) {
		printf("o");
		fflush(stdout);
	}

	if (0) {
error:
		ok = 0;
	}
	isl_set_free(set);
	isl_qpolynomial_free(cnt);
	isl_int_clear(v);
	isl_int_clear(n);
	isl_int_clear(d);
	isl_point_free(pnt);

	if (!ok)
		vpd->error = 1;

	if (vpd->options->continue_on_error)
		ok = 1;

	return (vpd->n >= 1 && ok) ? 0 : -1;
}

static int verify_isl(Polyhedron *P, Polyhedron *C,
		evalue *EP, const struct verify_options *options)
{
	struct verify_point_data vpd = { options };
	int i;
	isl_ctx *ctx = isl_ctx_alloc();
	isl_dim *dim;
	isl_set *set;
	isl_set *set_C;
	isl_int v;
	int r;

	dim = isl_dim_set_alloc(ctx, C->Dimension, P->Dimension - C->Dimension);
	for (i = 0; i < C->Dimension; ++i)
		dim = isl_dim_set_name(dim, isl_dim_param, i, options->params[i]);
	set = isl_set_new_from_polylib(P, isl_dim_copy(dim));
	dim = isl_dim_drop(dim, isl_dim_set, 0, P->Dimension - C->Dimension);
	set_C = isl_set_new_from_polylib(C, dim);
	set_C = isl_set_intersect(isl_set_copy(set), set_C);
	set_C = isl_set_remove(set_C, isl_dim_set, 0, P->Dimension - C->Dimension);

	set_C = set_bounds(set_C, options);

	isl_int_init(v);
	r = isl_set_count(set_C, &v);
	vpd.n = isl_int_cmp_si(v, 200) < 0 ? isl_int_get_si(v) : 200;
	isl_int_clear(v);

	if (!options->print_all) {
		vpd.s = vpd.n < 80 ? 1 : 1 + vpd.n/80;
		for (i = 0; i < vpd.n; i += vpd.s)
			printf(".");
		printf("\r");
		fflush(stdout);
	}

	vpd.set = set;
	vpd.pwqp = isl_pw_qpolynomial_from_evalue(isl_set_get_dim(set_C), EP);
	vpd.error = 0;
	if (r == 0)
		isl_set_foreach_point(set_C, verify_point, &vpd);
	if (vpd.error)
		r = -1;

	isl_pw_qpolynomial_free(vpd.pwqp);
	isl_set_free(set);
	isl_set_free(set_C);

	isl_ctx_free(ctx);

	if (!options->print_all)
		printf("\n");

	if (r < 0)
		fprintf(stderr, "Check failed !\n");

	return r;
}

static int verify(Polyhedron *P, Polyhedron *C, evalue *EP, skewed_gen_fun *gf,
		   arguments *options)
{
    Polyhedron *CS, *S;
    Vector *p;
    int result = 0;

    if (!options->series || options->function)
	return verify_isl(P, C, EP, &options->verify);

    CS = check_poly_context_scan(P, &C, C->Dimension, &options->verify);

    p = Vector_Alloc(P->Dimension+2);
    value_set_si(p->p[P->Dimension+1], 1);

    /* S = scanning list of polyhedra */
    S = Polyhedron_Scan(P, C, options->verify.barvinok->MaxRays);

    check_poly_init(C, &options->verify);

    /******* CHECK NOW *********/
    if (S) {
	if (!options->series || options->function) {
	    if (!check_poly_EP(S, CS, EP, 0, C->Dimension, 0, p->p,
				&options->verify))
		result = -1;
	} else {
	    if (!check_poly_gf(S, CS, gf, 0, C->Dimension, 0, p->p,
				&options->verify))
		result = -1;
	}
	Domain_Free(S);
    }

    if (result == -1)
	fprintf(stderr,"Check failed !\n");
    
    if (!options->verify.print_all)
	printf( "\n" );
  
    Vector_Free(p);
    if (CS) {
	Domain_Free(CS);
	Domain_Free(C);
    }

    return result;
}

/* frees M and Minv */
static void apply_transformation(Polyhedron **P, Polyhedron **C,
				 bool free_P, bool free_C,
				 Matrix *M, Matrix *Minv, Matrix **inv,
				 barvinok_options *options)
{
    Polyhedron *T;
    Matrix *M2;

    M2 = align_matrix(M, (*P)->Dimension + 1);
    T = *P;
    *P = Polyhedron_Preimage(*P, M2, options->MaxRays);
    if (free_P)
	Polyhedron_Free(T);
    Matrix_Free(M2);

    T = *C;
    *C = Polyhedron_Preimage(*C, M, options->MaxRays);
    if (free_C)
	Polyhedron_Free(T);

    Matrix_Free(M);

    if (*inv) {
	Matrix *T = *inv;
	*inv = Matrix_Alloc(Minv->NbRows, T->NbColumns);
	Matrix_Product(Minv, T, *inv);
	Matrix_Free(T);
	Matrix_Free(Minv);
    } else
	*inv = Minv;
}

/* Since we have "compressed" the parameters (in case there were
 * any equalities), the result is independent of the coordinates in the
 * coordinate subspace spanned by the lines.  We can therefore assume
 * these coordinates are zero and compute the inverse image of the map
 * from a lower dimensional space that adds zeros in the appropriate
 * places.
 */
static void remove_lines(Polyhedron *C, Matrix **M, Matrix **Minv)
{
    Matrix *L = Matrix_Alloc(C->Dimension+1, C->Dimension+1);
    for (int r = 0; r < C->NbBid; ++r)
	Vector_Copy(C->Ray[r]+1, L->p[r], C->Dimension);
    unimodular_complete(L, C->NbBid);
    assert(value_one_p(L->p[C->Dimension][C->Dimension]));
    assert(First_Non_Zero(L->p[C->Dimension], C->Dimension) == -1);
    Matrix_Transposition(L);
    assert(First_Non_Zero(L->p[C->Dimension], C->Dimension) == -1);

    *M = Matrix_Alloc(C->Dimension+1, C->Dimension-C->NbBid+1);
    for (int i = 0; i < C->Dimension+1; ++i)
	Vector_Copy(L->p[i]+C->NbBid, (*M)->p[i], C->Dimension-C->NbBid+1);

    Matrix *Linv = Matrix_Alloc(C->Dimension+1, C->Dimension+1);
    int ok = Matrix_Inverse(L, Linv);
    assert(ok);
    Matrix_Free(L);

    *Minv = Matrix_Alloc(C->Dimension-C->NbBid+1, C->Dimension+1);
    for (int i = C->NbBid; i < C->Dimension+1; ++i)
	Vector_AntiScale(Linv->p[i], (*Minv)->p[i-C->NbBid],
			 Linv->p[C->Dimension][C->Dimension], C->Dimension+1);
    Matrix_Free(Linv);
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

	if (C2->NbBid) {
	    Matrix *M, *Minv;
	    remove_lines(C2, &M, &Minv);
	    apply_transformation(&PT, &C2, PT != P, C2 != C, M, Minv, &inv,
				 options);
	}
    }
    POL_ENSURE_VERTICES(C2);
    if (!Polyhedron_has_revlex_positive_rays(C2, C2->Dimension)) {
	Polyhedron *T;
	Matrix *Constraints;
	Matrix *H, *Q, *U;
	Constraints = Matrix_Alloc(C2->NbConstraints, C2->Dimension+1);
	for (int i = 0; i < C2->NbConstraints; ++i)
	    Vector_Copy(C2->Constraint[i]+1, Constraints->p[i], C2->Dimension);
	left_hermite(Constraints, &H, &Q, &U);
	Matrix_Free(Constraints);
	/* flip rows of Q */
	for (int i = 0; i < C2->Dimension/2; ++i)
	    Vector_Exchange(Q->p[i], Q->p[C2->Dimension-1-i], C2->Dimension);
	Matrix_Free(H);
	Matrix_Free(U);
	Matrix *M = Matrix_Alloc(C2->Dimension+1, C2->Dimension+1);
	U = Matrix_Copy(Q);
	int ok = Matrix_Inverse(U, M);
	assert(ok);
	Matrix_Free(U);

	apply_transformation(&PT, &C2, PT != P, C2 != C, M, Q, &inv, options);
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
    skewed_gen_fun *gf = NULL;
    const char **param_name;
    int print_solution = 1;
    int result = 0;
    struct arguments options;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  		0 },
	{ &verify_argp,    	0,	"verification",		BV_GRP_LAST+1 },
	{ &convert_argp,    	0,	"output conversion",    BV_GRP_LAST+2 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    options.verify.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    M = Matrix_Read();
    assert(M);
    A = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    M = Matrix_Read();
    assert(M);
    C = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    assert(A->Dimension >= C->Dimension);
    param_name = Read_ParamNames(stdin, C->Dimension);

    if (options.verify.verify) {
	verify_options_set_range(&options.verify, A->Dimension);
	if (!bv_options->verbose)
	    print_solution = 0;
    }

    if (print_solution && bv_options->verbose) {
	Polyhedron_Print(stdout, P_VALUE_FMT, A);
	Polyhedron_Print(stdout, P_VALUE_FMT, C);
    }

    if (options.series) {
	gf = series(A, C, bv_options);
	if (print_solution) {
	    gf->print(cout, C->Dimension, param_name);
	    puts("");
	}
	if (options.function) {
	    EP = *gf;
	    if (print_solution)
		print_evalue(stdout, EP, param_name);
	}
    } else {
	EP = barvinok_enumerate_with_options(A, C, bv_options);
	assert(EP);
	if (evalue_convert(EP, &options.convert, bv_options->verbose,
			   C->Dimension, param_name))
	    print_solution = 0;
	if (options.size)
	    printf("\nSize: %d\n", evalue_size(EP));
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
    }

    if (options.verify.verify) {
	options.verify.params = param_name;
	result = verify(A, C, EP, gf, &options);
    }

    if (gf)
	delete gf;
    if (EP)
	evalue_free(EP);

    if (options.verify.barvinok->print_stats)
	barvinok_stats_print(options.verify.barvinok->stats, stdout);

    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    barvinok_options_free(bv_options);
    return result;
}
