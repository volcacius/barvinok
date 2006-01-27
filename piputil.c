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
    options->Simplify = 1;
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
