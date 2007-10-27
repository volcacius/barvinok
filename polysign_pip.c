#include <assert.h>
#include <piplib/piplibMP.h>
#include "polysign.h"

static PipQuast *solve_lp(int maximize, Matrix *C, Value *f, Value denom)
{
    int i;
    PipMatrix	*domain;
    PipOptions	*options;
    PipQuast   	*sol;

    domain = pip_matrix_alloc(C->NbRows+1, C->NbColumns+1);
    value_set_si(domain->p[0][0], 0);
    value_set_si(domain->p[0][1], -1);
    Vector_Copy(f, domain->p[0]+2, C->NbColumns-1);
    for (i = 0; i < C->NbRows; ++i) {
	value_assign(domain->p[i+1][0], C->p[i][0]);
	Vector_Copy(C->p[i]+1, domain->p[i+1]+2, C->NbColumns-1);
    }

    options = pip_options_init();
    options->Urs_unknowns = -1;
    options->Maximize = maximize;
    options->Nq = 0;
    sol = pip_solve(domain, NULL, -1, options);
    pip_matrix_free(domain);
    pip_options_free(options);

    return sol;
}

static enum lp_result constraints_affine_minmax(int maximize, Matrix *C,
					      Value *f, Value denom, Value *opt)
{
    enum lp_result res = lp_ok;
    PipQuast *sol = solve_lp(maximize, C, f, denom);

    if (!sol->list)
	res = lp_empty;
    else if (value_zero_p(sol->list->vector->the_deno[0]))
	res = lp_unbounded;
    else {
	if (maximize)
	    mpz_fdiv_q(*opt, sol->list->vector->the_vector[0],
			     sol->list->vector->the_deno[0]);
	else
	    mpz_cdiv_q(*opt, sol->list->vector->the_vector[0],
			     sol->list->vector->the_deno[0]);
    }
    pip_quast_free(sol);
    return res;
}

enum lp_result pip_constraints_opt(Matrix *C, Value *obj, Value denom,
				    enum lp_dir dir, Value *opt)
{
    int maximize = dir == lp_min ? 0 : 1;
    return constraints_affine_minmax(maximize, C, obj, denom, opt);
}

enum lp_result pip_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    enum lp_result res;
    Matrix M;

    if (emptyQ(D))
	return lp_empty;

    Polyhedron_Matrix_View(D, &M, D->NbConstraints);
    res = constraints_affine_minmax(0, &M, obj, denom, min);
    if (res != lp_ok)
	return res;
    res = constraints_affine_minmax(1, &M, obj, denom, max);
    return res;
}
