#include <assert.h>
#include <limits.h>
#include <sys/times.h>
#include <barvinok/barvinok.h>
#include "verify.h"
#include "argp.h"
#include "progname.h"

struct {
    int	    sign;
    int	    method;
    int	    flags;
} methods[] = {
{ BV_APPROX_SIGN_NONE,      BV_APPROX_NONE,	0 },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_BERNOULLI,	0 },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_DROP,	0 },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_VOLUME,	BV_VOL_LIFT },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_VOLUME,	BV_VOL_VERTEX },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_VOLUME,	BV_VOL_BARYCENTER },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_SCALE,	0 },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_SCALE,	BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST },
{ BV_APPROX_SIGN_APPROX,    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_DROP,	0 },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_VOLUME,	BV_VOL_LIFT },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_VOLUME,	BV_VOL_VERTEX },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_VOLUME,	BV_VOL_BARYCENTER },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	0 },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW  | BV_APPROX_SCALE_CHAMBER},
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW2 },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW2 | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW2 },
{ BV_APPROX_SIGN_LOWER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW2 | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_DROP,	0 },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_VOLUME,	BV_VOL_LIFT },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_VOLUME,	BV_VOL_VERTEX },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_VOLUME,	BV_VOL_BARYCENTER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	0 },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW2 },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_NARROW2 | BV_APPROX_SCALE_CHAMBER },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW2 },
{ BV_APPROX_SIGN_UPPER,	    BV_APPROX_SCALE,	BV_APPROX_SCALE_FAST | BV_APPROX_SCALE_NARROW2 | BV_APPROX_SCALE_CHAMBER },
};

#define nr_methods (sizeof(methods)/sizeof(*methods))

struct argp_option argp_options[] = {
    { "quiet",	    'q' },
    { 0 }
};

struct options {
    int quiet;
    struct verify_options    verify;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->verify.barvinok;
	state->child_inputs[1] = &options->verify;
	options->quiet = 0;
	break;
    case 'q':
	options->quiet = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct result_data {
    Value		    n;
    double		    RE_sum[nr_methods];

    clock_t		    ticks[nr_methods];
    size_t		    size[nr_methods];
};

void result_data_init(struct result_data *result)
{
    int i;
    for (i = 0; i < nr_methods; ++i) {
	result->RE_sum[i] = 0;
	result->ticks[i] = 0;
	result->size[i] = 0;
    }
    value_init(result->n);
}

void result_data_clear(struct result_data *result)
{
    int i;
    value_clear(result->n);
}

void result_data_print(struct result_data *result, int n)
{
    int i;

    fprintf(stderr, "%d", result->ticks[0]/n);
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %d", result->ticks[i]/n);
    fprintf(stderr, "\n");

    fprintf(stderr, "%d", result->size[0]/n);
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %d", result->size[i]/n);
    fprintf(stderr, "\n");

    fprintf(stderr, "%g\n", VALUE_TO_DOUBLE(result->n));
    fprintf(stderr, "%g", result->RE_sum[0]/VALUE_TO_DOUBLE(result->n));
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %g", result->RE_sum[i]/VALUE_TO_DOUBLE(result->n));
    fprintf(stderr, "\n");
}

struct test_approx_data {
    struct check_poly_data   cp;
    evalue	    	    **EP;
    struct result_data	    *result;
};

static void eval(const evalue *EP, Value *z, int sign, Value *v)
{
    evalue *res;

    res = evalue_eval(EP, z);
    if (sign == BV_APPROX_SIGN_LOWER)
	mpz_cdiv_q(*v, res->x.n, res->d);
    else if (sign == BV_APPROX_SIGN_UPPER)
	mpz_fdiv_q(*v, res->x.n, res->d);
    else if (sign == BV_APPROX_SIGN_APPROX)
	mpz_tdiv_q(*v, res->x.n, res->d);
    else {
	assert(value_one_p(res->d));
	value_assign(*v, res->x.n);
    }
    evalue_free(res);
}

static int test_approx(const struct check_poly_data *data, int nparam, Value *z,
		       const struct verify_options *options)
{
    struct test_approx_data* ta_data = (struct test_approx_data*) data;
    Value exact, approx;
    int i;

    value_init(exact);
    value_init(approx);

    eval(ta_data->EP[0], z, BV_APPROX_SIGN_NONE, &exact);

    /*
    value_print(stderr, VALUE_FMT, exact);
    */

    value_increment(ta_data->result->n, ta_data->result->n);
    for (i = 1; i < nr_methods; ++i) {
	double error;
	eval(ta_data->EP[i], z, methods[i].sign, &approx);
	/*
	fprintf(stderr, ", ");
	value_print(stderr, VALUE_FMT, approx);
	*/
	if (methods[i].sign == BV_APPROX_SIGN_LOWER)
	    assert(value_le(approx, exact));
	if (methods[i].sign == BV_APPROX_SIGN_UPPER)
	    assert(value_ge(approx, exact));
	value_subtract(approx, approx, exact);
	if (value_zero_p(exact))
	    error = abs(VALUE_TO_DOUBLE(approx));
	else
	    error = abs(VALUE_TO_DOUBLE(approx)) / VALUE_TO_DOUBLE(exact);
	ta_data->result->RE_sum[i] += error;
    }

    /*
    fprintf(stderr, "\n");
    */

    value_clear(exact);
    value_clear(approx);
    return 1;
}

static void test(Polyhedron *P, Polyhedron *C, evalue **EP,
		 struct result_data *result,
		 struct verify_options *options)
{
    Polyhedron *CS;
    Vector *p;
    unsigned nparam = C->Dimension;
    struct test_approx_data data;

    CS = check_poly_context_scan(P, &C, C->Dimension, options);

    p = Vector_Alloc(P->Dimension+2);
    value_set_si(p->p[P->Dimension+1], 1);

    check_poly_init(C, options);

    data.cp.z = p->p;
    data.cp.check = test_approx;
    data.EP = EP;
    data.result = result;
    check_poly(CS, &data.cp, nparam, 0, p->p+P->Dimension-nparam+1,
	       options);
    if (!options->print_all)
	printf("\n");

    Vector_Free(p);
    if (CS) {
	Domain_Free(CS);
	Domain_Free(C);
    }
}

void Matrix_File_Read_Input(FILE *in, Matrix *Mat)
{
  Value *p;
  int i,j,n;
  char *c, s[1024],str[1024];
  
  p = Mat->p_Init;
  for (i=0;i<Mat->NbRows;i++) {
    do {
      c = fgets(s, 1024, in);
      while(isspace(*c) && *c!='\n')
	++c;
    } while(c && (*c =='#' || *c== '\n'));
    
    if (!c) {
      errormsg1( "Matrix_Read", "baddim", "not enough rows" );
      break;
    }
    for (j=0;j<Mat->NbColumns;j++) {
      if(!c || *c=='\n' || *c=='#') {
	errormsg1("Matrix_Read", "baddim", "not enough columns");
	break;
      }
      if (sscanf(c,"%s%n",str,&n) == 0) {
	errormsg1( "Matrix_Read", "baddim", "not enough columns" );
	break;
      }
      value_read(*(p++),str);
      c += n;
    }
  }
} /* Matrix_Read_Input */

/* 
 * Read the contents of the matrix 'Mat' from standard input. 
 * A '#' in the first column is a comment line 
 */
Matrix *Matrix_File_Read(FILE *in)
{
  Matrix *Mat;
  unsigned NbRows, NbColumns;
  char s[1024];
  
  fgets(s, 1024, in);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d %d", &NbRows, &NbColumns)<2))
    fgets(s, 1024, in);
  Mat = Matrix_Alloc(NbRows,NbColumns);
  if(!Mat) {
    errormsg1("Matrix_Read", "outofmem", "out of memory space");
    return(NULL);
  }
  Matrix_File_Read_Input(in, Mat);
  return Mat;
} /* Matrix_Read */

void handle(FILE *in, struct result_data *result, struct verify_options *options)
{
    int i;
    Polyhedron *A, *C;
    Matrix *M;
    const char **param_name;
    evalue *EP[nr_methods];

    M = Matrix_File_Read(in);
    A = Constraints2Polyhedron(M, options->barvinok->MaxRays);
    Matrix_Free(M);
    M = Matrix_File_Read(in);
    C = Constraints2Polyhedron(M, options->barvinok->MaxRays);
    Matrix_Free(M);
    param_name = Read_ParamNames(in, C->Dimension);

    for (i = 0; i < nr_methods; ++i) {
	struct tms st_cpu;
	struct tms en_cpu;
	options->barvinok->polynomial_approximation = methods[i].sign;
	options->barvinok->approximation_method = methods[i].method;
	if (methods[i].method == BV_APPROX_SCALE)
	    options->barvinok->scale_flags = methods[i].flags;
	else if (methods[i].method == BV_APPROX_VOLUME)
	    options->barvinok->volume_triangulate = methods[i].flags;

	times(&st_cpu);
	EP[i] = barvinok_enumerate_with_options(A, C, options->barvinok);
	times(&en_cpu);
	result->ticks[i] = en_cpu.tms_utime - st_cpu.tms_utime;
	/*
	print_evalue(stdout, EP[i], param_name);
	*/
    }
    for (i = 0; i < nr_methods; ++i)
	result->size[i] = evalue_size(EP[i])/4;
    test(A, C, EP, result, options);
    for (i = 0; i < nr_methods; ++i)
	evalue_free(EP[i]);

    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
}

int main(int argc, char **argv)
{
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    char path[PATH_MAX+1];
    struct result_data all_result;
    int n = 0;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  		0 },
	{ &verify_argp,    	0,	"verification",		BV_GRP_LAST+1 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct options options;

    options.verify.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    if (options.verify.M == INT_MIN)
	options.verify.M = 100;
    if (options.verify.m == INT_MAX)
	options.verify.m = -100;

    result_data_init(&all_result);

    while (fgets(path, sizeof(path), stdin)) {
	struct result_data result;
	FILE *in;
	int i;

	++n;
	result_data_init(&result);
	fprintf(stderr, "%s", path);
	*strchr(path, '\n') = '\0';
	in = fopen(path, "r");
	assert(in);
	handle(in, &result, &options.verify);
	fclose(in);

	if (!options.quiet)
	    result_data_print(&result, 1);

	value_addto(all_result.n, all_result.n, result.n);
	for (i = 0; i < nr_methods; ++i) {
	    all_result.RE_sum[i] += result.RE_sum[i];
	    all_result.ticks[i] += result.ticks[i];
	    all_result.size[i] += result.size[i];
	}

	result_data_clear(&result);

	if (!options.quiet) {
	    fprintf(stderr, "average\n");
	    result_data_print(&all_result, n);
	}
    }

    result_data_clear(&all_result);

    barvinok_options_free(bv_options);

    return 0;
}
