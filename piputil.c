#include <polylib/polylibgmp.h>
#include <stdio.h>  /* needed for piplib ! */
#include <piplib/piplibMP.h>
#include <assert.h>
#include "piputil.h"
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

static PipMatrix *poly2pip(Polyhedron *P, int pos, int n, int nparam)
{
    PipMatrix *	    matrix;
    int i, j;
    unsigned int extra = P->Dimension - pos - n - nparam;

    matrix = pip_matrix_alloc(P->NbConstraints, P->Dimension+2);
    for (i = 0; i < P->NbConstraints; ++i) {
	value_assign(matrix->p[i][0], P->Constraint[i][0]);
	for (j = 0; j < n; ++j)
	    value_assign(matrix->p[i][1+j], P->Constraint[i][1+pos+j]);
	for (j = 0; j < pos; ++j)
	    value_assign(matrix->p[i][1+n+j], P->Constraint[i][1+j]);
	for (j = 0; j < extra; ++j)
	    value_assign(matrix->p[i][1+n+pos+j], P->Constraint[i][1+n+pos+j]);
	for (j = 1+pos+n+extra; j < P->Dimension+2; ++j)
	    value_assign(matrix->p[i][j], P->Constraint[i][j]);
    }

    return matrix;
}

static int max_new(PipQuast *q, int max, int d, int *maxd)
{
    PipNewparm *p;

    for (p = q->newparm; p; p = p->next)
	if (p->rank > max)
	    max = p->rank;
    if (q->condition) {
	if (++d > *maxd)
	    *maxd = d;
	max = max_new(q->next_else, max, d, maxd);
	max = max_new(q->next_then, max, d, maxd);
    }
    return max;
}

/*
 * v: total number of variables in original problem
 * pos: position of minimized variables
 * n: dimension of space in which minimum was taken
 *    i.e., number of minimized vars in original problem
 * e: number of extra existential variables
 * d: total dimension
 */
struct scan_data {
    Matrix *M;
    int	    v;
    int	    pos;
    int	    n;
    int	    e;
    int	    d;
    void    (*f)(struct scan_data *data, PipQuast *q, int i);
};

static int vectorpos2cpos(struct scan_data *data, int j)
{
    int np = data->d - data->v - data->e;
    int l;

    if (j < data->pos)
	l = 1+j;
    else if (j < data->v - data->n)
	l = 1+data->n+j;
    else if (j < data->v-data->n+np)
	l = 1+data->n+data->e+j;
    else
	l = 1+data->n+j-np;

    return l;
}

/*
 * i: next row
 *
 * Dimensions in the resulting domain are:
 *	v n e (d-v-n-e)
 */
static void scan_quast_r(struct scan_data *data, PipQuast *q, int i)
{
    PipNewparm *p;

    for (p = q->newparm; p; p = p->next) {
	int j, l;
	PipVector *vec = p->vector;

	value_set_si(data->M->p[i][0], 1);
	Vector_Set(data->M->p[i]+1, 0, data->d);

	for (j = 0; j < vec->nb_elements-1; ++j) {
	    l = vectorpos2cpos(data, j);
	    value_assign(data->M->p[i][l], vec->the_vector[j]);
	}
	l = vectorpos2cpos(data, p->rank);
	value_oppose(data->M->p[i][l], p->deno);
	value_assign(data->M->p[i][1+data->d], vec->the_vector[j]);
	++i;

	value_set_si(data->M->p[i][0], 1);
	Vector_Set(data->M->p[i]+1, 0, data->d);

	for (j = 0; j < vec->nb_elements-1; ++j) {
	    l = vectorpos2cpos(data, j);
	    value_oppose(data->M->p[i][l], vec->the_vector[j]);
	}
	l = vectorpos2cpos(data, p->rank);
	value_assign(data->M->p[i][l], p->deno);
	value_oppose(data->M->p[i][1+data->d], vec->the_vector[j]);
	value_addto(data->M->p[i][1+data->d], data->M->p[i][1+data->d], p->deno);
	value_decrement(data->M->p[i][1+data->d], data->M->p[i][1+data->d]);
	++i;
    }

    if (q->condition) {
	int j;
	value_set_si(data->M->p[i][0], 1);
	Vector_Set(data->M->p[i]+1, 0, data->d);

	for (j = 0; j < q->condition->nb_elements-1; ++j) {
	    int l = vectorpos2cpos(data, j);
	    value_assign(data->M->p[i][l], q->condition->the_vector[j]);
	}
	value_assign(data->M->p[i][1+data->d], q->condition->the_vector[j]);
	scan_quast_r(data, q->next_then, i+1);

	for (j = 0; j < q->condition->nb_elements-1; ++j) {
	    int l = vectorpos2cpos(data, j);
	    value_oppose(data->M->p[i][l], q->condition->the_vector[j]);
	}
	value_oppose(data->M->p[i][1+data->d], q->condition->the_vector[j]);
	value_decrement(data->M->p[i][1+data->d], data->M->p[i][1+data->d]);
	scan_quast_r(data, q->next_else, i+1);

	return;
    }

    /* if data->n is zero, we are only interested in the domains
     * where a solution exists and not in the actual solution
     */
    if (q->list && data->n) {
	PipList *l;
	int j, k;
	for (j = 0, l = q->list; l; ++j, l = l->next) {
	    Vector_Set(data->M->p[i], 0, data->d+1);
	    value_set_si(data->M->p[i][1+data->pos+j], -1);

	    for (k = 0; k < l->vector->nb_elements-1; ++k) {
		int ll = vectorpos2cpos(data, k);
		value_assign(data->M->p[i][ll], l->vector->the_vector[k]);
	    }
	    value_assign(data->M->p[i][1+data->d], l->vector->the_vector[k]);

	    ++i;
	}
    }
    data->f(data, q, i);
}

static void scan_quast(struct scan_data *data, PipQuast *q)
{
    int		nparam, nexist, d, i, rows;

    nparam = data->d - data->n;

    d = 0;
    nexist = max_new(q, nparam-1, 0, &d) - nparam+1;
    rows = nparam + 2 * nexist + d + data->n;

    /* nparam now refers to the number of parameters in the original polyhedron */
    nparam -= data->v - data->n;

    data->e = nexist;
    data->d = data->v + nexist + nparam;

    data->M = Matrix_Alloc(rows, 1+data->d+1);
    /* All parameters are/should be positive */
    for (i = 0; i < data->pos; ++i) {
	value_set_si(data->M->p[i][0], 1);
	value_set_si(data->M->p[i][1+i], 1);
    }
    for (i = data->pos+data->n; i < data->v; ++i) {
	value_set_si(data->M->p[i-data->n][0], 1);
	value_set_si(data->M->p[i-data->n][1+i], 1);
    }
    for (i = 0; i < nparam; ++i) {
	value_set_si(data->M->p[data->v-data->n+i][0], 1);
	value_set_si(data->M->p[data->v-data->n+i][1+data->v+nexist+i], 1);
    }

    scan_quast_r(data, q, data->d - data->e - data->n);
    Matrix_Free(data->M);
}

struct quast2poly_data {
    struct scan_data	scan;
    Polyhedron*		D;
};

static void add_quast_base(struct scan_data *data, PipQuast *q, int i)
{
    struct quast2poly_data *qdata = (struct quast2poly_data *)data;
    if (q->list) {
	Matrix *C;
	Polyhedron *P;
	C = Matrix_Alloc(i, 1+data->d+1);
	Vector_Copy(data->M->p[0], C->p[0], i * (1+data->d+1));
	P = Constraints2Polyhedron(C, MAXRAYS);
	Matrix_Free(C);
	P->next = qdata->D;
	qdata->D = P;
    }
}

/*
 * nvar: total number of variables, including the minimized variables,
 *	 but excluding the parameters
 * nparam: the number of parameters in the PIP problem, excluding the "big parameter"
 *	 for maximization problems, i.e., the
 *	 total number of variables minus the minimized variables
 * pos: position of the first minimized variable
 * n: number of minimized variables
 */
static Polyhedron *quast2poly(PipQuast *q, int nvar, int nparam, int pos, int n)
{
    struct quast2poly_data	data;

    data.scan.v = nvar;
    data.scan.pos = pos;
    data.scan.n = n;
    data.scan.d = n+nparam;

    data.D = 0;
    data.scan.f = add_quast_base;

    scan_quast(&data.scan, q);

    return data.D;
}

Polyhedron *pip_lexminmax(Polyhedron *P, int pos, int n, int nparam, int max)
{
    PipOptions	*options;
    PipMatrix *	domain;
    unsigned int nvar = P->Dimension - nparam;
    PipMatrix   *context = pip_matrix_alloc(0, P->Dimension - n + 2);
    PipQuast	*sol;
    Polyhedron  *min;

    POL_ENSURE_INEQUALITIES(P);

    domain = poly2pip(P, pos, n, nparam);
#ifdef PIP_DEBUG
    pip_matrix_print(stderr, domain);
    pip_matrix_print(stderr, context);
#endif

    options = pip_options_init();
    options->Maximize = max;
    options->Urs_unknowns = -1;
    sol = pip_solve(domain, context, -1, options);

#ifdef PIP_DEBUG
    pip_quast_print(stderr, sol, 0);
#endif

    min = quast2poly(sol, nvar, P->Dimension - n, pos, n);

    pip_quast_free(sol);
    pip_matrix_free(context);
    pip_matrix_free(domain);
    pip_options_free(options);

    return min;
}

Polyhedron *pip_projectout(Polyhedron *P, int pos, int n, int nparam)
{
    PipOptions	*options;
    PipMatrix *	domain;
    unsigned int nvar = P->Dimension - nparam;
    PipMatrix   *context = pip_matrix_alloc(0, P->Dimension - n + 2);
    PipQuast	*sol;
    Polyhedron  *min;

    POL_ENSURE_INEQUALITIES(P);

    domain = poly2pip(P, pos, n, nparam);
#ifdef PIP_DEBUG
    pip_matrix_print(stderr, domain);
    pip_matrix_print(stderr, context);
#endif

    options = pip_options_init();
    sol = pip_solve(domain, context, -1, options);

#ifdef PIP_DEBUG
    pip_quast_print(stderr, sol, 0);
#endif

    min = quast2poly(sol, nvar - n, P->Dimension - n, 0, 0);

    pip_quast_free(sol);
    pip_matrix_free(context);
    pip_matrix_free(domain);
    pip_options_free(options);

    return min;
}

/* The variables in the mapping are 
 *	previous_user next_user extra parameters
 * previous_user vars are constrained using equalities in M
 *
 * The variables in the remaining domains are
 *	first_user extra parameters
 */
struct inputdep_data {
    struct scan_data	 scan;
    Polyhedron*		 D;	    /* The remaining domains */
    Polyhedron*		 M;	    /* The mapping */
    PipOptions		*options;
    PipMatrix   	*context;
    PipMatrix 		*domain;
    int			 offset;
    int		     	 depth;
};

static int quast_empty(PipQuast *q)
{
    return !q->newparm && !q->condition && !q->list;
}

static void inputdep_at(struct inputdep_data *data)
{
    PipQuast	*sol;

    value_set_si(data->domain->p[data->offset+data->depth][0], 1);
    value_set_si(data->domain->p[data->offset+data->depth]
				[data->domain->NbColumns-1], -1);

    sol = pip_solve(data->domain, data->context, -1, data->options);

    scan_quast(&data->scan, sol);
    pip_quast_free(sol);
}

static void inputdep_base(struct scan_data *data, PipQuast *q, int row)
{
    struct inputdep_data *idata = (struct inputdep_data *)data;
    if (q->list) {
	Matrix *C;
	Polyhedron *P;
	C = Matrix_Alloc(row, 1+data->d+1);
	Vector_Copy(data->M->p[0], C->p[0], row * (1+data->d+1));
	P = Constraints2Polyhedron(C, MAXRAYS);
	Matrix_Free(C);
	P->next = idata->M;
	idata->M = P;
    } else if (idata->depth > 0) {
	int i, j;
	struct inputdep_data newdata;
	/* Ignore the positivity constraints */
	int offset = data->d - data->e - data->n;
	int extra = row - offset;

	if (data->e != 0)
	    newdata.context = pip_matrix_alloc(0, data->d - data->n + 2);
	else
	    newdata.context = idata->context;

	newdata.domain = pip_matrix_alloc(idata->offset + idata->depth + extra,
					  data->d + 2);

	for (i = 0; i < idata->offset; ++i) {
	    for (j = 0; j <= data->v; ++j)
		value_assign(newdata.domain->p[i][j], idata->domain->p[i][j]);
	    for ( ; j <= data->d - data->e + 1; ++j)
		value_assign(newdata.domain->p[i][j+data->e], 
			     idata->domain->p[i][j]);
	}
	for (i = 0; i < extra; ++i)
	    for (j = 0; j <= data->d+1; ++j)
		value_assign(newdata.domain->p[idata->offset + i][j],
			     data->M->p[offset + i][j]);
	for (i = 0; i < idata->depth; ++i) {
	    for (j = 0; j <= data->v; ++j)
		value_assign(newdata.domain->p[idata->offset+extra + i][j], 
			     idata->domain->p[idata->offset + i][j]);
	    for ( ; j <= data->d - data->e + 1; ++j)
		value_assign(newdata.domain->p[idata->offset+extra + i][j+data->e], 
			     idata->domain->p[idata->offset + i][j]);
	}

	newdata.D = idata->D;
	newdata.M = idata->M;
	newdata.options = idata->options;
	newdata.offset = idata->offset + extra;
	newdata.depth = idata->depth-1;

	newdata.scan.v = data->v + data->e;
	newdata.scan.pos = data->pos;
	newdata.scan.n = data->n;
	newdata.scan.d = data->d;
	newdata.scan.f = inputdep_base;

	inputdep_at(&newdata);

	idata->D = newdata.D;
	idata->M = newdata.M;

	pip_matrix_free(newdata.domain);
	if (data->e != 0)
	    pip_matrix_free(newdata.context);
    } else {
	int i, j;
	PipMatrix *	domain;
	PipQuast	*sol;
	PipOptions	*options;
	PipMatrix   *context;
	Polyhedron	*P;
	int offset = data->d - data->e - data->n;
	int extra = row - offset;

	if (data->e != 0)
	    context = pip_matrix_alloc(0, data->d - data->n + 2);
	else
	    context = idata->context;

	domain = pip_matrix_alloc(idata->offset + extra + data->n,
				  data->d + 2);
	/* It's ok (though unneeded) here to also copy the constraints
	 * that impose that the two iterations access the same data
	 * since we will be imposing the (stronger) constraints that
	 * the two iterations are the same below.
	 */
	for (i = 0; i < idata->offset; ++i) {
	    for (j = 0; j <= data->v; ++j)
		value_assign(domain->p[i][j], idata->domain->p[i][j]);
	    for ( ; j <= data->d - data->e + 1; ++j)
		value_assign(domain->p[i][j+data->e], 
			     idata->domain->p[i][j]);
	}
	for (i = 0; i < extra; ++i)
	    for (j = 0; j <= data->d+1; ++j)
		value_assign(domain->p[idata->offset + i][j],
			     data->M->p[offset + i][j]);
	for (i = 0; i < data->n; ++i) {
	    value_set_si(domain->p[idata->offset + extra + i][1 + i], 1);
	    value_set_si(domain->p[idata->offset + extra + i][1 + data->n + i], -1);
	}
	options = pip_options_init();
	sol = pip_solve(domain, context, -1, options);

	P = quast2poly(sol, data->v + data->e, data->d - data->n, data->pos, data->n);
	if (P) {
	    Matrix *C;
	    assert(P->next == NULL);
	    /* The constraints contain data->n equalities equating the
	     * iterators of the two domains
	     */
	    C = Matrix_Alloc(P->NbConstraints - data->n, P->Dimension - data->n + 2);
	    for (i = j = 0; i < P->NbConstraints; ++i) {
		if (First_Non_Zero(&P->Constraint[i][1], data->n) != -1)
		    continue;
		value_assign(C->p[j][0], P->Constraint[i][0]);
		Vector_Copy(P->Constraint[i] + 1+data->n, C->p[j] + 1, 
			    C->NbColumns-1);
		++j;
	    }
	    assert(i - j == data->n);
	    Polyhedron_Free(P);
	    P = Constraints2Polyhedron(C, MAXRAYS);
	    Matrix_Free(C);
	    P->next = idata->D;
	    idata->D = P;
	}

	pip_matrix_free(domain);
	if (data->e != 0)
	    pip_matrix_free(context);
	pip_options_free(options);
    }
}

struct pip_dep_result pip_inputdep(Polyhedron *D, int dim, Matrix *M)
{
    struct pip_dep_result res;
    struct inputdep_data data;
    int		 nparam = D->Dimension - dim;
    int i, j, k;

    PipQuast	*sol;

    data.context = pip_matrix_alloc(0, D->Dimension + 2);
    data.domain = pip_matrix_alloc(2*D->NbConstraints + M->NbRows + dim, 
				   2 + dim + D->Dimension);
    for (i = 0; i < 2; ++i) {
	int d = i*D->NbConstraints;
	for (j = 0; j < D->NbConstraints; ++j) {
	    value_assign(data.domain->p[d+j][0], D->Constraint[j][0]);
	    for (k = 1; k <= dim; ++k)
		value_assign(data.domain->p[d+j][i*dim+k], D->Constraint[j][k]);
	    for ( ; k < D->Dimension+2; ++k)
		value_assign(data.domain->p[d+j][dim+k], D->Constraint[j][k]);
	}
    }
    data.offset = 2*D->NbConstraints;
    for (j = 0; j < M->NbRows; ++j)
	for (k = 0; k < dim; ++k) {
	    value_oppose(data.domain->p[data.offset+j][1+k], M->p[j][k]);
	    value_assign(data.domain->p[data.offset+j][1+dim+k], M->p[j][k]);
	}
    data.offset += M->NbRows;
    for (k = 0; k < dim; ++k) {
	value_set_si(data.domain->p[data.offset+k][1+k], -1);
	value_set_si(data.domain->p[data.offset+k][1+dim+k], 1);
    }

    data.options = pip_options_init();
    data.options->Maximize = 1;
    data.options->Simplify = 1;

    data.scan.v = 2 * dim;
    data.scan.pos = 0;
    data.scan.n = dim;
    data.scan.d = dim + D->Dimension;

    data.D = NULL;
    data.M = NULL;
    data.scan.f = inputdep_base;
    data.depth = dim-1;

    inputdep_at(&data);

    pip_matrix_free(data.context);
    pip_matrix_free(data.domain);
    pip_options_free(data.options);

    res.D = data.D;
    res.M = data.M;

    return res;
}
