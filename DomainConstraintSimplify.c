#include <polylib/polylibgmp.h>

/*
 * Replaces constraint a x >= c by x >= ceil(c/a)
 * where "a" is a common factor in the coefficients
 * old is the constraint; v points to an initialized
 * value that this procedure can use.
 * Return non-zero if something changed.
 * Result is placed in new.
 */
int ConstraintSimplify(Value *old, Value *new, int len, Value* v)
{
    Vector_Gcd(old+1, len - 2, v);

    if (value_one_p(*v))
	return 0;

    Vector_AntiScale(old+1, new+1, *v, len-2);
    mpz_fdiv_q(new[len-1], old[len-1], *v);
    return 1;
}

static Polyhedron *p_simplify_constraints(Polyhedron *P, Vector *row,
					  Value *g, unsigned MaxRays)
{
    Polyhedron *T, *R = P;
    int len = P->Dimension+2;
    int r;

    /* Also look at equalities.
     * If an equality can be "simplified" then there
     * are no integer solutions anyway and the following loop
     * will add a conflicting constraint
     */
    for (r = 0; r < R->NbConstraints; ++r) {
	if (ConstraintSimplify(R->Constraint[r], row->p, len, g)) {
	    T = R;
	    R = AddConstraints(row->p, 1, R, MaxRays);
	    if (T != P)
		Polyhedron_Free(T);
	    r = -1;
	}
    }
    if (R != P)
	Polyhedron_Free(P);
    return R;
}

/*
 * Replaces constraint a x >= c by x >= ceil(c/a)
 * where "a" is a common factor in the coefficients
 * Destroys P and returns a newly allocated Polyhedron
 * or just returns P in case no changes were made
 */
Polyhedron *DomainConstraintSimplify(Polyhedron *P, unsigned MaxRays)
{
    Polyhedron **prev;
    int len = P->Dimension+2;
    Vector *row = Vector_Alloc(len);
    Value g;
    Polyhedron *R = P, *N;
    value_set_si(row->p[0], 1);
    value_init(g);

    for (prev = &R; P; P = N) {
	Polyhedron *T;
	N = P->next;
	T = p_simplify_constraints(P, row, &g, MaxRays);

	if (emptyQ(T) && prev != &R) {
	    Polyhedron_Free(T);
	    *prev = NULL;
	    continue;
	}

	if (T != P)
	    T->next = N;
	*prev = T;
	prev = &T->next;
    }

    if (R->next && emptyQ(R)) {
	N = R->next;
	Polyhedron_Free(R);
	R = N;
    }

    value_clear(g);
    Vector_Free(row);
    return R;
}
