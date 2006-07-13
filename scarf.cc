#include <vector>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "scarf.h"

using std::vector;

static Matrix *extract_matrix(Polyhedron *P, unsigned dim)
{
    Matrix *A;
    int n_col;

    n_col = 0;
    for (int i = 0; i < P->NbConstraints; ++i)
	if (value_notzero_p(P->Constraint[i][1+dim]) ||
	    value_notzero_p(P->Constraint[i][1+dim+1]))
	    ++n_col;

    assert(n_col == 3);

    A = Matrix_Alloc(2, n_col+2);
    n_col = 0;
    for (int i = 0; i < P->NbConstraints; ++i) {
	if (value_zero_p(P->Constraint[i][1+dim]) &&
	    value_zero_p(P->Constraint[i][1+dim+1]))
	    continue;
	value_assign(A->p[0][n_col], P->Constraint[i][1+dim]);
	value_assign(A->p[1][n_col], P->Constraint[i][1+dim+1]);
	++n_col;
    }
    value_set_si(A->p[0][n_col], 1);
    value_set_si(A->p[1][n_col+1], 1);

    return A;
}

static int lex_sign(Value *v, int len)
{
    int first;

    first = First_Non_Zero(v, len);
    return first == -1 ? 0 : value_sign(v[first]);
}

static Matrix *normalize_matrix(Matrix *A, int pos[4], int n)
{
    Matrix *T, *B;
    Value tmp, factor;
    int type;

    value_init(tmp);
    value_init(factor);

    T = Matrix_Alloc(2, 2);
    Extended_Euclid(A->p[0][pos[0]], A->p[1][pos[0]], 
		    &T->p[0][0], &T->p[0][1], &tmp);
    value_division(T->p[1][0], A->p[1][pos[0]], tmp);
    value_division(T->p[1][1], A->p[0][pos[0]], tmp);
    value_oppose(T->p[0][0], T->p[0][0]);
    value_oppose(T->p[0][1], T->p[0][1]);
    value_oppose(T->p[1][0], T->p[1][0]);

    B = Matrix_Alloc(2, A->NbColumns);
    Matrix_Product(T, A, B);
    Matrix_Free(T);

    if (lex_sign(B->p[1], B->NbColumns) > 0) {
	value_set_si(tmp, -1);
	Vector_Scale(B->p[1], B->p[1], tmp, B->NbColumns);
    }

    assert(n == 3);
    if (n == 3) {
	if (value_notzero_p(B->p[1][pos[1]]) &&
	    value_notzero_p(B->p[1][pos[2]])) {
	    assert(value_neg_p(B->p[1][pos[1]]));
	    assert(value_pos_p(B->p[1][pos[2]]));
	    assert(value_pos_p(B->p[0][pos[1]]));
	    assert(value_neg_p(B->p[0][pos[2]]));
	    Vector_Exchange(B->p[0], B->p[1], B->NbColumns);
	    type = 2;
	} else {
	    assert(0);
	}
    }

    assert(type == 2);

    if (type == 2) {
	for (;;) {
	    value_oppose(tmp, B->p[0][pos[1]]);
	    value_pdivision(factor, tmp, B->p[1][pos[1]]);
	    value_oppose(tmp, B->p[1][pos[2]]);
	    value_pdivision(tmp, tmp, B->p[0][pos[2]]);
	    if (value_zero_p(factor) && value_zero_p(tmp))
		break;
	    assert(value_zero_p(factor) || value_zero_p(tmp));
	    if (value_pos_p(factor)) {
		value_set_si(tmp, 1);
		Vector_Combine(B->p[0], B->p[1], B->p[0], tmp, factor, B->NbColumns);
		if (value_zero_p(B->p[0][pos[1]])) {
		    /* We will deal with this later */
		    assert(lex_sign(B->p[0], B->NbColumns) < 0);
		}
	    } else {
		value_set_si(factor, 1);
		Vector_Combine(B->p[0], B->p[1], B->p[1], tmp, factor, B->NbColumns);
		if (value_zero_p(B->p[1][pos[2]])) {
		    /* We will deal with this later */
		    assert(lex_sign(B->p[1], B->NbColumns) < 0);
		}
	    }
	}
    }

    value_clear(tmp);
    value_clear(factor);

    return B;
}

struct simplex {
    Matrix * M;

    simplex(int d) {
	M = Matrix_Alloc(d, 2);
    }
    void transform(Matrix *T);
    Polyhedron *shrunk_polyhedron(Polyhedron *P, int dim, Matrix *A,
				    unsigned MaxRays);
};

static bool lex_smaller(Value *v1, Value *v2, int n)
{
    for (int i = 0; i < n; ++i)
	if (value_lt(v1[i], v2[i]))
	    return true;
	else if (value_gt(v1[i], v2[i]))
	    return false;
    return false;
}

void simplex::transform(Matrix *T)
{
    Matrix *M2 = M;
    M = Matrix_Alloc(M2->NbRows, M2->NbColumns);
    Matrix_Product(M2, T, M);
    Matrix_Free(M2);

    int lexmin = 0;
    for (int i = 1; i < M->NbRows; ++i)
	if (lex_smaller(M->p[i], M->p[lexmin], 2))
	    lexmin = i;
    if (lex_sign(M->p[lexmin], 2) < 0) {
	Value tmp;
	value_init(tmp);
	value_set_si(tmp, -1);
	Vector_Scale(M->p[lexmin], M->p[lexmin], tmp, 2);
	value_set_si(tmp, 1);
	for (int i = 0; i < M->NbRows; ++i) {
	    if (i == lexmin)
		continue;
	    Vector_Combine(M->p[lexmin], M->p[i], M->p[i], tmp, tmp, 2);
	}
	value_clear(tmp);
    }
}

Polyhedron *simplex::shrunk_polyhedron(Polyhedron *P, int dim, Matrix *A,
					unsigned MaxRays)
{
    Matrix *Constraints, *offset;
    Polyhedron *Q;
    Value min;
    value_init(min);

    offset = Matrix_Alloc(M->NbRows, A->NbColumns);
    Matrix_Product(M, A, offset);

    Constraints = Polyhedron2Constraints(P);

    for (int i = 0, j = 0; i < Constraints->NbRows; ++i) {
	if (value_zero_p(Constraints->p[i][1+dim]) &&
	    value_zero_p(Constraints->p[i][1+dim+1]))
	    continue;
	value_set_si(min, 0);
	for (int k = 0; k < offset->NbRows; ++k)
	    if (value_lt(offset->p[k][j], min))
		value_assign(min, offset->p[k][j]);
	value_addto(Constraints->p[i][1+P->Dimension],
		    Constraints->p[i][1+P->Dimension], min);
	++j;
    }
    Q = Constraints2Polyhedron(Constraints, MaxRays);

    Matrix_Free(offset);
    Matrix_Free(Constraints);
    value_clear(min);

    return Q;
}

struct scarf_complex {
    vector<simplex> simplices;
    void add(Matrix *B, int pos[4], int n);
    void print(FILE *out);
    ~scarf_complex() {
	for (int i = 0; i < simplices.size(); ++i)
	    Matrix_Free(simplices[i].M);
    }
};

void scarf_complex::add(Matrix *B, int pos[4], int n)
{
    Matrix *T;

    assert(n == 3);

    T = Matrix_Alloc(2, 2);
    Vector_Copy(B->p[0]+B->NbColumns-2, T->p[0], 2);
    Vector_Copy(B->p[1]+B->NbColumns-2, T->p[1], 2);

    if (n == 3) {
	simplex s1(1);
	value_set_si(s1.M->p[0][0], 0);
	value_set_si(s1.M->p[0][1], 1);
	s1.transform(T);
	simplices.push_back(s1);

	simplex s2(1);
	value_set_si(s2.M->p[0][0], 1);
	value_set_si(s2.M->p[0][1], 1);
	s2.transform(T);
	simplices.push_back(s2);

	simplex s3(1);
	value_set_si(s3.M->p[0][0], 1);
	value_set_si(s3.M->p[0][1], 0);
	s3.transform(T);
	simplices.push_back(s3);

	simplex s4(2);
	value_set_si(s4.M->p[0][0], 0);
	value_set_si(s4.M->p[0][1], 1);
	value_set_si(s4.M->p[1][0], 1);
	value_set_si(s4.M->p[1][1], 1);
	s4.transform(T);
	simplices.push_back(s4);

	simplex s5(2);
	value_set_si(s5.M->p[0][0], 1);
	value_set_si(s5.M->p[0][1], 0);
	value_set_si(s5.M->p[1][0], 1);
	value_set_si(s5.M->p[1][1], 1);
	s5.transform(T);
	simplices.push_back(s5);
    }

    Matrix_Free(T);
}

void scarf_complex::print(FILE *out)
{
    for (int i = 0; i < simplices.size(); ++i) {
	Matrix_Print(out, P_VALUE_FMT, simplices[i].M);
    }
}

struct scarf_collector {
    virtual void add(Polyhedron *P, int sign, Polyhedron *C, unsigned MaxRays) = 0;
};

static void scarf(Polyhedron *P, unsigned exist, unsigned nparam, unsigned MaxRays,
		  scarf_collector& col)
{
    Matrix *A, *B;
    int dim = P->Dimension - exist - nparam;
    assert(exist == 2);
    int pos[4];
    Polyhedron *U;
    gen_fun *gf;

    A = extract_matrix(P, dim);

    for (int i = 0; i < 3; ++i)
	pos[i] = i;
    B = normalize_matrix(A, pos, 3);

    scarf_complex scarf;
    scarf.add(B, pos, 3);

    U = Universe_Polyhedron(nparam);
    col.add(P, 0, U, MaxRays);
    for (int i = 0; i < scarf.simplices.size(); ++i) {
	Polyhedron *Q;
	int sign = (scarf.simplices[i].M->NbRows % 2) ? -1 : 1;
	Q = scarf.simplices[i].shrunk_polyhedron(P, dim, A, MaxRays);
	col.add(Q, sign, U, MaxRays);
	Polyhedron_Free(Q);
    }
    Polyhedron_Free(U);

    Matrix_Free(B);

    Matrix_Free(A);
}

struct scarf_collector_gf : public scarf_collector {
    ZZ cn, cd;
    gen_fun *gf;

    scarf_collector_gf() {
	cd = 1;
    }
    virtual void add(Polyhedron *P, int sign, Polyhedron *C, unsigned MaxRays);
};

void scarf_collector_gf::add(Polyhedron *P, int sign, Polyhedron *C, 
			     unsigned MaxRays)
{
    if (!sign)
	gf = barvinok_series(P, C, MaxRays);
    else {
	gen_fun *gf2;
	cn = sign;
	gf2 = barvinok_series(P, C, MaxRays);
	gf->add(cn, cd, gf2);
	delete gf2;
    }
}

gen_fun *barvinok_enumerate_scarf_series(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    scarf_collector_gf scgf;
    scarf(P, exist, nparam, MaxRays, scgf);
    return scgf.gf;
}

struct scarf_collector_ev : public scarf_collector {
    evalue mone;
    evalue *EP;

    scarf_collector_ev() {
	value_init(mone.d);
	evalue_set_si(&mone, -1, 1);
    }
    ~scarf_collector_ev() {
	free_evalue_refs(&mone); 
    }
    virtual void add(Polyhedron *P, int sign, Polyhedron *C, unsigned MaxRays);
};

void scarf_collector_ev::add(Polyhedron *P, int sign, Polyhedron *C, 
			     unsigned MaxRays)
{
    if (!sign)
	EP = barvinok_enumerate_ev(P, C, MaxRays);
    else {
	evalue *E2;
	E2 = barvinok_enumerate_ev(P, C, MaxRays);
	if (sign < 0)
	    emul(&mone, E2);
	eadd(E2, EP);
	free_evalue_refs(E2);
	free(E2);
    }
}

evalue *barvinok_enumerate_scarf(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays)
{
    scarf_collector_ev scev;
    scarf(P, exist, nparam, MaxRays, scev);
    return scev.EP;
}
