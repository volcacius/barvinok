#include <assert.h>
#include <piplib/piplibMP.h>
#include <barvinok/basis_reduction.h>

struct pip_lp {
    Polyhedron *P;
    Value      *obj;
    Matrix     *eq;
    int	        neq;
    PipQuast   *sol;
};

static struct pip_lp *init_lp(Polyhedron *P);
static void set_lp_obj(struct pip_lp *lp, Value *row, int dim);
static int solve_lp(struct pip_lp *lp);
static void get_obj_val(struct pip_lp* lp, mpq_t *F);
static void delete_lp(struct pip_lp *lp);
static int add_lp_row(struct pip_lp *lp, Value *row, int dim);
static void get_alpha(struct pip_lp* lp, int row, mpq_t *alpha);
static void del_lp_row(struct pip_lp *lp);

#define GBR_LP			    	    struct pip_lp
#define GBR_type		    	    mpq_t
#define GBR_init(v)		    	    mpq_init(v)
#define GBR_clear(v)		    	    mpq_clear(v)
#define GBR_set(a,b)			    mpq_set(a,b)
#define GBR_set_ui(a,b)			    mpq_set_ui(a,b,1)
#define GBR_mul(a,b,c)			    mpq_mul(a,b,c)
#define GBR_lt(a,b)			    (mpq_cmp(a,b) < 0)
#define GBR_floor(a,b)			    mpz_fdiv_q(a,mpq_numref(b),mpq_denref(b))
#define GBR_ceil(a,b)			    mpz_cdiv_q(a,mpq_numref(b),mpq_denref(b))
#define GBR_lp_init(P)		    	    init_lp(P)
#define GBR_lp_set_obj(lp, obj, dim)	    set_lp_obj(lp, obj, dim)
#define GBR_lp_solve(lp)		    solve_lp(lp)
#define GBR_lp_get_obj_val(lp, F)	    get_obj_val(lp, F)
#define GBR_lp_delete(lp)		    delete_lp(lp)
#define GBR_lp_next_row(lp)		    lp->neq
#define GBR_lp_add_row(lp, row, dim)	    add_lp_row(lp, row, dim)
#define GBR_lp_get_alpha(lp, row, alpha)    get_alpha(lp, row, alpha)
#define GBR_lp_del_row(lp)		    del_lp_row(lp);
#define Polyhedron_Reduced_Basis    	    pip_dual_Polyhedron_Reduced_Basis
#include "basis_reduction_templ.c"

#define ALLOC(type) (type*)malloc(sizeof(type))

static struct pip_lp *init_lp(Polyhedron *P)
{
    struct pip_lp *lp = ALLOC(struct pip_lp);
    lp->P = P;
    lp->eq = Matrix_Alloc(P->Dimension, P->Dimension);
    lp->neq = 0;
    lp->sol = NULL;

    return lp;
}

static void set_lp_obj(struct pip_lp *lp, Value *row, int dim)
{
    assert(lp->P->Dimension == dim);
    lp->obj = row;
}

static int solve_lp(struct pip_lp *lp)
{
    PipOptions	*options;
    PipMatrix	*domain;
    unsigned	 dim = lp->P->Dimension;
    int		 rows = 1 + 2*dim + 2*lp->P->NbConstraints;
    int		 cols = 1 + 1 + lp->neq + 2*lp->P->NbConstraints + 1;
    int		 i, j;

    if (lp->sol)
	pip_quast_free(lp->sol);

    assert(lp->P->NbEq == 0);
    domain = pip_matrix_alloc(rows, cols);
    for (i = 0; i < lp->neq; ++i)
	for (j = 0; j < dim; ++j) {
	    value_assign(domain->p[j][1+1+i], lp->eq->p[i][j]);
	    value_oppose(domain->p[dim+j][1+1+i], lp->eq->p[i][j]);
	}
    value_set_si(domain->p[2*dim][1], 1);
    for (j = 0; j < dim; ++j) {
	value_oppose(domain->p[j][cols-1], lp->obj[j]);
	value_assign(domain->p[dim+j][cols-1], lp->obj[j]);
    }
    for (i = 0; i < lp->P->NbConstraints; ++i) {
	for (j = 0; j < dim; ++j) {
	    value_assign(domain->p[j][1+1+lp->neq+i], lp->P->Constraint[i][1+j]);
	    value_assign(domain->p[dim+j][1+1+lp->neq+lp->P->NbConstraints+i],
			 lp->P->Constraint[i][1+j]);
	}
	value_assign(domain->p[2*dim][1+1+lp->neq+i], lp->P->Constraint[i][1+dim]);
	value_assign(domain->p[2*dim][1+1+lp->neq+lp->P->NbConstraints+i],
		     lp->P->Constraint[i][1+dim]);
	value_set_si(domain->p[2*dim+1+i][0], 1);
	value_set_si(domain->p[2*dim+1+i][1+1+lp->neq+i], -1);
	value_set_si(domain->p[2*dim+1+lp->P->NbConstraints+i][0], 1);
	value_set_si(domain->p[2*dim+1+lp->P->NbConstraints+i]
				[1+1+lp->neq+lp->P->NbConstraints+i], -1);
    }

    options = pip_options_init();
    options->Urs_unknowns = -1;
    options->Nq = 0;
    lp->sol = pip_solve(domain, NULL, -1, options);
    pip_matrix_free(domain);

    pip_options_free(options);

    return value_zero_p(lp->sol->list->vector->the_deno[0]);
}

static void get_obj_val(struct pip_lp* lp, mpq_t *F)
{
    value_assign(mpq_numref(*F), lp->sol->list->vector->the_vector[0]);
    value_assign(mpq_denref(*F), lp->sol->list->vector->the_deno[0]);
}

static void delete_lp(struct pip_lp *lp)
{
    if (lp->sol)
	pip_quast_free(lp->sol);
    Matrix_Free(lp->eq);
    free(lp);
}

static int add_lp_row(struct pip_lp *lp, Value *row, int dim)
{
    assert(lp->P->Dimension == dim);
    Vector_Copy(row, lp->eq->p[lp->neq], dim);
    return lp->neq++;
}

static void get_alpha(struct pip_lp* lp, int row, mpq_t *alpha)
{
    PipList *list;
    int i = 0;

    for (i = 0, list = lp->sol->list->next; i < row; ++i, list = list->next)
	;

    if (value_zero_p(list->vector->the_deno[0])) {
	value_set_si(mpq_numref(*alpha), 0);
	value_set_si(mpq_denref(*alpha), 1);
    } else {
	value_oppose(mpq_numref(*alpha), list->vector->the_vector[0]);
	value_assign(mpq_denref(*alpha), list->vector->the_deno[0]);
    }
}

static void del_lp_row(struct pip_lp *lp)
{
    assert(lp->neq > 0);
    lp->neq--;
}
