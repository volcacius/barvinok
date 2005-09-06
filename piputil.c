#include <polylib/polylibgmp.h>
#include <stdio.h>  /* needed for piplib ! */
#include <piplib/piplibMP.h>
#include <assert.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

static PipMatrix *poly2pip(Polyhedron *P, int exist, int nparam)
{
    PipMatrix *	    matrix;
    int i, j;
    unsigned int nvar = P->Dimension - exist - nparam;

    matrix = pip_matrix_alloc(P->NbConstraints, P->Dimension+2);
    for (i = 0; i < P->NbConstraints; ++i) {
	value_assign(matrix->p[i][0], P->Constraint[i][0]);
	for (j = 0; j < exist; ++j)
	    value_assign(matrix->p[i][1+j], P->Constraint[i][1+nvar+j]);
	for (j = 0; j < nvar; ++j)
	    value_assign(matrix->p[i][1+exist+j], P->Constraint[i][1+j]);
	for (j = 1+exist+nvar; j < P->Dimension+2; ++j)
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

static int vectorpos2cpos(int v, int n, int e, int d, int j)
{
    int np = d-v-n-e;
    int l;

    if (j < v)
	l = 1+j;
    else if (j < np+v)
	l = 1+v+n+e+j-v;
    else
	l = 1+v+n+j-v-np;

    return l;
}

/*
 * i: next row
 * n: dimension of space in which minimum was taken
 *    i.e., number of existential vars in original problem
 * v: number of variables in original problem
 * e: number of extra existential variables
 * d: total dimension
 *
 * Dimensions in the resulting domain are:
 *	v n e (d-v-n-e)
 */
static void add_quast(Polyhedron**D, Matrix* M, PipQuast *q, 
		      int i, int v, int n, int e, int d)
{
    PipNewparm *p;

    for (p = q->newparm; p; p = p->next) {
	int j, l;
	PipVector *vec = p->vector;

	value_set_si(M->p[i][0], 1);
	Vector_Set(M->p[i]+1, 0, d);

	for (j = 0; j < vec->nb_elements-1; ++j) {
	    l = vectorpos2cpos(v, n, e, d, j);
	    value_assign(M->p[i][l], vec->the_vector[j]);
	}
	l = vectorpos2cpos(v, n, e, d, p->rank);
	value_oppose(M->p[i][l], p->deno);
	value_assign(M->p[i][1+d], vec->the_vector[j]);
	++i;

	value_set_si(M->p[i][0], 1);
	Vector_Set(M->p[i]+1, 0, d);

	for (j = 0; j < vec->nb_elements-1; ++j) {
	    l = vectorpos2cpos(v, n, e, d, j);
	    value_oppose(M->p[i][l], vec->the_vector[j]);
	}
	l = vectorpos2cpos(v, n, e, d, p->rank);
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
	    int l = vectorpos2cpos(v, n, e, d, j);
	    value_assign(M->p[i][l], q->condition->the_vector[j]);
	}
	value_assign(M->p[i][1+d], q->condition->the_vector[j]);
	add_quast(D, M, q->next_then, i+1, v, n, e, d);

	for (j = 0; j < q->condition->nb_elements-1; ++j) {
	    int l = vectorpos2cpos(v, n, e, d, j);
	    value_oppose(M->p[i][l], q->condition->the_vector[j]);
	}
	value_oppose(M->p[i][1+d], q->condition->the_vector[j]);
	value_decrement(M->p[i][1+d], M->p[i][1+d]);
	add_quast(D, M, q->next_else, i+1, v, n, e, d);
    } else if (q->list) {
	PipList *l;
	Matrix *C;
	Polyhedron *P;
	int j, k;
	for (j = 0, l = q->list; l; ++j, l = l->next) {
	    Vector_Set(M->p[i], 0, d+1);
	    value_set_si(M->p[i][1+v+j], -1);

	    for (k = 0; k < l->vector->nb_elements-1; ++k) {
		int ll = vectorpos2cpos(v, n, e, d, k);
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

static Polyhedron *quast2poly(PipQuast *q, int nvar, int nset)
{
    int			nparam, nexist;
    PipList*		l;
    int			i, d, rows;
    Matrix*		M;
    Polyhedron*		D;

    assert(nset != 0);	/* required ? */
    nparam = q->newparm ? q->newparm->rank :
	     q->condition ? q->condition->nb_elements-1 :
			    q->list->vector->nb_elements-1;
    d = 0;
    nexist = max_new(q, nparam-1, 0, &d) - nparam+1;
    rows = nparam + 2 * nexist + d + nset;
    M = Matrix_Alloc(rows, 1+nset+nexist+nparam+1);
    /* All parameters are/should be positive */
    for (i = 0; i < nvar; ++i) {
	value_set_si(M->p[i][0], 1);
	value_set_si(M->p[i][1+i], 1);
    }
    for ( ; i < nparam; ++i) {
	value_set_si(M->p[i][0], 1);
	value_set_si(M->p[i][1+nvar+nset+nexist+i-nvar], 1);
    }
    D = 0;
    add_quast(&D, M, q, nparam, nvar, nset, nexist, nset+nexist+nparam);
    Matrix_Free(M);
    return D;
}

Polyhedron *pip_lexmin(Polyhedron *P, int exist, int nparam)
{
    PipOptions	*options;
    PipMatrix *	domain;
    unsigned int nvar = P->Dimension - exist - nparam;
    PipMatrix   *context = pip_matrix_alloc(0, nvar+nparam+2);
    PipQuast	*sol;
    Polyhedron  *min;

    POL_ENSURE_INEQUALITIES(P);

    domain = poly2pip(P, exist, nparam);
    /*
    pip_matrix_print(stderr, domain);
    pip_matrix_print(stderr, context);
    */

    options = pip_options_init();
    sol = pip_solve(domain, context, -1, options);

    /*
    pip_quast_print(stderr, sol, 0);
    */

    min = quast2poly(sol, nvar, exist);

    pip_quast_free(sol);
    pip_matrix_free(context);
    pip_matrix_free(domain);
    pip_options_free(options);

    return min;
}
