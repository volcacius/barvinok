#include <sstream>
#include "fdstream.h"
#include <barvinok/util.h>
#include "edomain.h"
#include "evalue_util.h"

using std::vector;
using std::endl;
using std::ostream;

static void print_term(ostream& os, Value v, int pos, int dim,
		        char **names, int *first)
{
    if (value_zero_p(v)) {
	if (first && *first && pos >= dim)
	    os << "0";
	return;
    }

    if (first) {
	if (!*first && value_pos_p(v))
	    os << "+";
	*first = 0;
    }
    if (pos < dim) {
	if (value_mone_p(v)) {
	    os << "-";
	} else if (!value_one_p(v))
	    os << VALUE_TO_INT(v);
	os << names[pos];
    } else
	os << VALUE_TO_INT(v);
}

void EDomain_floor::print(ostream& os, char **p) const
{
    int first = 1;
    os << "[";
    os << "(";
    for (int i = 0; i < v->Size-2; ++i)
	print_term(os, v->p[1+i], i, v->Size-2, p, &first);
    os << ")";
    os << "/";
    print_term(os, v->p[0], v->Size-2, v->Size-2, p, NULL);
    os << "]";
}

void EDomain::print_constraints(ostream& os, char **p,
				barvinok_options *options) const
{
    Value tmp;
    value_init(tmp);

    Matrix *M = Matrix_Alloc(2*floors.size(), 1+D->Dimension+1);
    value_set_si(tmp, -1);
    for (int i = 0; i < floors.size(); ++i) {
	value_set_si(M->p[2*i][0], 1);
	Vector_Copy(floors[i]->v->p+1, M->p[2*i]+1, dimension());
	value_assign(M->p[2*i][1+D->Dimension], floors[i]->v->p[1+dimension()]);
	value_oppose(M->p[2*i][1+dimension()+i], floors[i]->v->p[0]);

	Vector_Scale(M->p[2*i]+1, M->p[2*i+1]+1, tmp, D->Dimension+1);
	value_addto(M->p[2*i+1][1+D->Dimension], M->p[2*i+1][1+D->Dimension],
		    M->p[2*i+1][1+dimension()+i]);
	value_decrement(M->p[2*i+1][1+D->Dimension], M->p[2*i+1][1+D->Dimension]);
	value_set_si(M->p[2*i+1][0], 1);
    }
    Polyhedron *E = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    Polyhedron *SD = DomainSimplify(D, E, options->MaxRays);
    Polyhedron_Free(E);

    char **names = p;
    unsigned dim = dimension();
    if (dim < SD->Dimension) {
	names = new char * [SD->Dimension];
	int i;
	for (i = 0; i < dim; ++i)
	    names[i] = p[i];
	for ( ; i < SD->Dimension; ++i) {
	    std::ostringstream strm;
	    floors[i-dim]->print(strm, p);
	    names[i] = strdup(strm.str().c_str());
	}
    }

    for (int i = 0; i < SD->NbConstraints; ++i) {
	int first = 1;
	int v = First_Non_Zero(SD->Constraint[i]+1, SD->Dimension);
	if (v == -1)
	    continue;
	if (i)
	    os << " && ";
	if (value_pos_p(SD->Constraint[i][v+1])) {
	    print_term(os, SD->Constraint[i][v+1], v, SD->Dimension,
		       names, NULL);
	    if (value_zero_p(SD->Constraint[i][0]))
		os << " = ";
	    else
		os << " >= ";
	    for (int j = v+1; j <= SD->Dimension; ++j) {
		value_oppose(tmp, SD->Constraint[i][1+j]);
		print_term(os, tmp, j, SD->Dimension,
			   names, &first);
	    }
	} else {
	    value_oppose(tmp, SD->Constraint[i][1+v]);
	    print_term(os, tmp, v, SD->Dimension,
		       names, NULL);
	    if (value_zero_p(SD->Constraint[i][0]))
		os << " = ";
	    else
		os << " <= ";
	    for (int j = v+1; j <= SD->Dimension; ++j)
		print_term(os, SD->Constraint[i][1+j], j, SD->Dimension,
			   names, &first);
	}
    }
    value_clear(tmp);
    Domain_Free(SD);

    if (dim < D->Dimension) {
	for (int i = dim; i < D->Dimension; ++i)
	    free(names[i]);
	delete [] names;
    }
}

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

EDomain_floor::EDomain_floor(const evalue *f, int dim)
{
    e = new evalue;
    value_init(e->d);
    evalue_copy(e, f);
    v = Vector_Alloc(1+dim+1);
    value_set_si(v->p[0], 1);
    evalue2constraint_r(NULL, e, v->p, v->Size);
    refcount = 1;
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

	    EDomain_floor *new_floor;
	    new_floor = new EDomain_floor(&e->x.p->arr[0], dimension());

	    /* constraints for the new floor */
	    int row = D->NbConstraints+2*fract;
	    Vector_Copy(new_floor->v->p+1, M->p[row]+1, dimension());
	    value_assign(M->p[row][cols-1], new_floor->v->p[1+dimension()]);
	    value_oppose(M->p[row][1+D->Dimension+fract], new_floor->v->p[0]);
	    value_set_si(M->p[row][0], 1);
	    assert(value_eq(M->p[row][cols-1], new_floor->v->p[1+dimension()]));
	    assert(Vector_Equal(new_floor->v->p+1, M->p[row]+1, dimension()));

	    Vector_Scale(M->p[row]+1, M->p[row+1]+1, mone.x.n, cols-1);
	    value_set_si(M->p[row+1][0], 1);
	    value_addto(M->p[row+1][cols-1], M->p[row+1][cols-1],
			M->p[row+1][1+D->Dimension+fract]);
	    value_decrement(M->p[row+1][cols-1], M->p[row+1][cols-1]);

	    new_floors.push_back(new_floor);

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
