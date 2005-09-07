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

static int vectorpos2cpos(int v, int pos, int n, int e, int d, int j)
{
    int np = d-v-e;
    int l;

    if (j < pos)
	l = 1+j;
    else if (j < v-n)
	l = 1+n+j;
    else if (j < v-n+np)
	l = 1+n+e+j;
    else
	l = 1+n+j-np;

    return l;
}

/*
 * i: next row
 * v: total number of variables in original problem
 * pos: position of minimized variables
 * n: dimension of space in which minimum was taken
 *    i.e., number of minimized vars in original problem
 * e: number of extra existential variables
 * d: total dimension
 *
 * Dimensions in the resulting domain are:
 *	v n e (d-v-n-e)
 */
static void add_quast(Polyhedron**D, Matrix* M, PipQuast *q, 
		      int i, int v, int pos, int n, int e, int d)
{
    PipNewparm *p;

    for (p = q->newparm; p; p = p->next) {
	int j, l;
	PipVector *vec = p->vector;

	value_set_si(M->p[i][0], 1);
	Vector_Set(M->p[i]+1, 0, d);

	for (j = 0; j < vec->nb_elements-1; ++j) {
	    l = vectorpos2cpos(v, pos, n, e, d, j);
	    value_assign(M->p[i][l], vec->the_vector[j]);
	}
	l = vectorpos2cpos(v, pos, n, e, d, p->rank);
	value_oppose(M->p[i][l], p->deno);
	value_assign(M->p[i][1+d], vec->the_vector[j]);
	++i;

	value_set_si(M->p[i][0], 1);
	Vector_Set(M->p[i]+1, 0, d);

	for (j = 0; j < vec->nb_elements-1; ++j) {
	    l = vectorpos2cpos(v, pos, n, e, d, j);
	    value_oppose(M->p[i][l], vec->the_vector[j]);
	}
	l = vectorpos2cpos(v, pos, n, e, d, p->rank);
	value_assign(M->p[i][l], p->deno);
	value_oppose(M->p[i][1+d], vec->the_vector[j]);
	value_addto(M->p[i][1+d], M->p[i][1+d], p->deno);
	value_decrement(M->p[i][1+d], M->p[i][1+d]);
	++i;
    }

    if (q->condition) {
	int j;
	value_set_si(M->p[i][0], 1);
	Vector_Set(M->p[i]+1, 0, d);

	for (j = 0; j < q->condition->nb_elements-1; ++j) {
	    int l = vectorpos2cpos(v, pos, n, e, d, j);
	    value_assign(M->p[i][l], q->condition->the_vector[j]);
	}
	value_assign(M->p[i][1+d], q->condition->the_vector[j]);
	add_quast(D, M, q->next_then, i+1, v, pos, n, e, d);

	for (j = 0; j < q->condition->nb_elements-1; ++j) {
	    int l = vectorpos2cpos(v, pos, n, e, d, j);
	    value_oppose(M->p[i][l], q->condition->the_vector[j]);
	}
	value_oppose(M->p[i][1+d], q->condition->the_vector[j]);
	value_decrement(M->p[i][1+d], M->p[i][1+d]);
	add_quast(D, M, q->next_else, i+1, v, pos, n, e, d);
    } else if (q->list) {
	PipList *l;
	Matrix *C;
	Polyhedron *P;
	int j, k;
	for (j = 0, l = q->list; l; ++j, l = l->next) {
	    Vector_Set(M->p[i], 0, d+1);
	    value_set_si(M->p[i][1+pos+j], -1);

	    for (k = 0; k < l->vector->nb_elements-1; ++k) {
		int ll = vectorpos2cpos(v, pos, n, e, d, k);
		value_assign(M->p[i][ll], l->vector->the_vector[k]);
	    }
	    value_assign(M->p[i][1+d], l->vector->the_vector[k]);

	    ++i;
	}
	C = Matrix_Alloc(i, 1+d+1);
	Vector_Copy(M->p[0], C->p[0], i * (1+d+1));
	P = Constraints2Polyhedron(C, MAXRAYS);
	Matrix_Free(C);
	P->next = *D;
	*D = P;
    }
}

/*
 * nvar: total number of variables, including the minimized variables,
 *	 but excluding the parameters
 * pos: position of the first minimized variable
 * n: number of minimized variables
 */
static Polyhedron *quast2poly(PipQuast *q, int nvar, int pos, int n)
{
    int			nparam;
    int			nexist;
    PipList*		l;
    int			i, d, rows;
    Matrix*		M;
    Polyhedron*		D;

    assert(n != 0);	/* required ? */
    nparam = q->newparm ? q->newparm->rank :
	     q->condition ? q->condition->nb_elements-1 :
			    q->list->vector->nb_elements-1;
    d = 0;
    nexist = max_new(q, nparam-1, 0, &d) - nparam+1;
    rows = nparam + 2 * nexist + d + n;

    /* nparam now refers to the number of parameters in the original polyhedron */
    nparam -= nvar - n;
    M = Matrix_Alloc(rows, 1+nvar+nexist+nparam+1);
    /* All parameters are/should be positive */
    for (i = 0; i < pos; ++i) {
	value_set_si(M->p[i][0], 1);
	value_set_si(M->p[i][1+i], 1);
    }
    for (i = pos+n; i < nvar; ++i) {
	value_set_si(M->p[i-n][0], 1);
	value_set_si(M->p[i-n][1+i], 1);
    }
    for (i = 0; i < nparam; ++i) {
	value_set_si(M->p[nvar-n+i][0], 1);
	value_set_si(M->p[nvar-n+i][1+nvar+nexist+i], 1);
    }
    D = 0;
    add_quast(&D, M, q, nparam, nvar, pos, n, nexist, nvar+nexist+nparam);
    Matrix_Free(M);
    return D;
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
    options->Max = max;
    sol = pip_solve(domain, context, -1, options);

#ifdef PIP_DEBUG
    pip_quast_print(stderr, sol, 0);
#endif

    min = quast2poly(sol, nvar, pos, n);

    pip_quast_free(sol);
    pip_matrix_free(context);
    pip_matrix_free(domain);
    pip_options_free(options);

    return min;
}
