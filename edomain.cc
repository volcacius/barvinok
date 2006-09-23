#include "fdstream.h"
#include <barvinok/util.h>
#include "edomain.h"
#include "evalue_util.h"

using std::vector;
using std::endl;

void EDomain::print(FILE *out, char **p)
{
    fdostream os(dup(fileno(out)));
    for (int i = 0; i < floors.size(); ++i) {
	os << "floor " << i << ": [";
	evalue_print(os, floors[i]->e, p);
	os << "]" << endl;
    }
    Polyhedron_Print(out, P_VALUE_FMT, D);
}

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 : 
	  p->type == flooring ? 1 : 0;
}

static void add_coeff(Value *cons, int len, evalue *coeff, int pos)
{
    Value tmp;

    assert(value_notzero_p(coeff->d));

    value_init(tmp);

    value_lcm(cons[0], coeff->d, &tmp);
    value_division(tmp, tmp, cons[0]);
    Vector_Scale(cons, cons, tmp, len);
    value_division(tmp, cons[0], coeff->d);
    value_addmul(cons[pos], tmp, coeff->x.n);

    value_clear(tmp);
}

static int evalue2constraint_r(EDomain *D, evalue *E, Value *cons, int len);

static void add_fract(evalue *e, Value *cons, int len, int pos)
{
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    /* contribution of alpha * fract(X) is 
     *      alpha * X ...
     */
    assert(e->x.p->size == 3);
    evalue argument;
    value_init(argument.d);
    evalue_copy(&argument, &e->x.p->arr[0]);
    emul(&e->x.p->arr[2], &argument);
    evalue2constraint_r(NULL, &argument, cons, len);
    free_evalue_refs(&argument);

    /*	    - alpha * floor(X) */
    emul(&mone, &e->x.p->arr[2]);
    add_coeff(cons, len, &e->x.p->arr[2], pos);
    emul(&mone, &e->x.p->arr[2]);

    free_evalue_refs(&mone); 
}

static int evalue2constraint_r(EDomain *D, evalue *E, Value *cons, int len)
{
    int r = 0;
    if (value_zero_p(E->d) && E->x.p->type == fractional) {
	int i;
	assert(E->x.p->size == 3);
	r = evalue2constraint_r(D, &E->x.p->arr[1], cons, len);
	assert(value_notzero_p(E->x.p->arr[2].d));
	if (D && (i = D->find_floor(&E->x.p->arr[0])) >= 0) {
	    add_fract(E, cons, len, 1+D->D->Dimension-D->floors.size()+i);
	} else {
	    if (value_pos_p(E->x.p->arr[2].x.n)) {
		evalue coeff;
		value_init(coeff.d);
		value_init(coeff.x.n);
		value_set_si(coeff.d, 1);
		evalue_denom(&E->x.p->arr[0], &coeff.d);
		value_decrement(coeff.x.n, coeff.d);
		emul(&E->x.p->arr[2], &coeff);
		add_coeff(cons, len, &coeff, len-1);
		free_evalue_refs(&coeff);
	    }
	    r = 1;
	}
    } else if (value_zero_p(E->d)) {
	assert(E->x.p->type == polynomial);
	assert(E->x.p->size == 2);
	r = evalue2constraint_r(D, &E->x.p->arr[0], cons, len) || r;
	add_coeff(cons, len, &E->x.p->arr[1], E->x.p->pos);
    } else {
	add_coeff(cons, len, E, len-1);
    }
    return r;
}

int evalue2constraint(EDomain *D, evalue *E, Value *cons, int len)
{
    Vector_Set(cons, 0, len);
    value_set_si(cons[0], 1);
    return evalue2constraint_r(D, E, cons, len);
}

Matrix *EDomain::add_ge_constraint(evalue *constraint,
				   vector<EDomain_floor *>& new_floors) const
{
    evalue mone;
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
    int fract = 0;
    for (evalue *e = constraint; value_zero_p(e->d); 
	 e = &e->x.p->arr[type_offset(e->x.p)]) {
	int i;
	if (e->x.p->type != fractional)
	    continue;
	if (find_floor(&e->x.p->arr[0]) == -1)
	    ++fract;
    }

    int rows = D->NbConstraints+2*fract+1;
    int cols = 2+D->Dimension+fract;
    Matrix *M = Matrix_Alloc(rows, cols);
    for (int i = 0; i < D->NbConstraints; ++i) {
	Vector_Copy(D->Constraint[i], M->p[i], 1+D->Dimension);
	value_assign(M->p[i][1+D->Dimension+fract], 
		     D->Constraint[i][1+D->Dimension]);
    }
    value_set_si(M->p[rows-1][0], 1);
    fract = 0;
    evalue *e;
    for (e = constraint; value_zero_p(e->d); e = &e->x.p->arr[type_offset(e->x.p)]) {
	if (e->x.p->type == fractional) {
	    int i, pos;

	    i = find_floor(&e->x.p->arr[0]);
	    if (i >= 0)
		pos = D->Dimension-floors.size()+i;
	    else
		pos = D->Dimension+fract;

	    add_fract(e, M->p[rows-1], cols, 1+pos);

	    if (pos < D->Dimension)
		continue;

	    /* constraints for the new floor */
	    int row = D->NbConstraints+2*fract;
	    value_set_si(M->p[row][0], 1);
	    evalue2constraint_r(NULL, &e->x.p->arr[0], M->p[row], cols);
	    value_oppose(M->p[row][1+D->Dimension+fract], M->p[row][0]);
	    value_set_si(M->p[row][0], 1);

	    Vector_Scale(M->p[row]+1, M->p[row+1]+1, mone.x.n, cols-1);
	    value_set_si(M->p[row+1][0], 1);
	    value_addto(M->p[row+1][cols-1], M->p[row+1][cols-1],
			M->p[row+1][1+D->Dimension+fract]);
	    value_decrement(M->p[row+1][cols-1], M->p[row+1][cols-1]);

	    new_floors.push_back(new EDomain_floor(&e->x.p->arr[0]));

	    ++fract;
	} else {
	    assert(e->x.p->type == polynomial);
	    assert(e->x.p->size == 2);
	    add_coeff(M->p[rows-1], cols, &e->x.p->arr[1], e->x.p->pos);
	}
    }
    add_coeff(M->p[rows-1], cols, e, cols-1);
    value_set_si(M->p[rows-1][0], 1);
    free_evalue_refs(&mone); 
    return M;
}
