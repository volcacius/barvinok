#include <iostream>
#include <barvinok/genfun.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include <barvinok/basis_reduction.h>
#include "config.h"

using std::cout;
using std::cerr;
using std::endl;

static Polyhedron *uncone(Polyhedron *C, unsigned MaxRays)
{
    Matrix *Constraints;
    Polyhedron *P;
    int i;

    Constraints = Matrix_Alloc(C->NbConstraints, 1+(C->Dimension-1)+1);
    for (i = 0; i < C->NbConstraints; ++i) {
	/* positivity constraints */
	if (First_Non_Zero(C->Constraint[i]+1, C->Dimension) == -1) {
	    assert(i == C->NbConstraints-1);
	    Constraints->NbRows = i;
	    break;
	}
	assert(value_zero_p(C->Constraint[i][1+C->Dimension]));
	Vector_Copy(C->Constraint[i], Constraints->p[i], 1+C->Dimension);
    }
    P = Constraints2Polyhedron(Constraints, MaxRays);
    Matrix_Free(Constraints);

    return P;
}

void static scan(gen_fun *S)
{
    Value c;
    Vector *params;
    value_init(c);
    params = Vector_Alloc(2);
    for (int j = 19; j <= 39; ++j) {
	value_set_si(params->p[1], j);
	for (int i = 8; i <= 17; ++i) {
	    value_set_si(params->p[0], i);
	    S->coefficient(params->p, &c);
	    fprintf(stderr, "%d %d: ", i, j);
	    value_print(stderr, VALUE_FMT, c);
	    fprintf(stderr, "\n");
	}
    }
    value_clear(c);
    Vector_Free(params);
}

bool is_indicator(gen_fun *gf, barvinok_options *options)
{
    gen_fun *test;
    QQ mone(-1, 1);

    test = gf->Hadamard_product(gf, options);
    test->add(mone, gf);

    test->print(std::cerr, 0, NULL);
    cerr << endl;
}

int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *P, *C, *B, *Pold;
    Vector *small;
    Matrix *T, *basis;
    gen_fun *S, *S_shift, *hp, *S_divide, *frob;
    vec_ZZ up;
    Value c;
    bool ok;
    QQ mone(-1, 1);
    Vector *coins;
    barvinok_options *options = barvinok_options_new_with_defaults();

    coins = Vector_Read();
    assert(coins);
    assert(coins->Size == 4);

    M = Matrix_Alloc(5, 7);
    Vector_Copy(coins->p, M->p[0]+1, 4);
    Vector_Free(coins);
    value_set_si(M->p[0][1+4], -1);
    for (int i = 0; i < 4; ++i) {
	value_set_si(M->p[1+i][0], 1);
	value_set_si(M->p[1+i][1+i], 1);
    }
    C = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    C = remove_equalities_p(C, 4, NULL);
    assert(C->NbEq == 0);
    Polyhedron_Print(stderr, P_VALUE_FMT, C);

    B = uncone(C, options->MaxRays);

    basis = Polyhedron_Reduced_Basis(B, options);
    small = Vector_Alloc(B->Dimension + 2);
    Vector_Copy(basis->p[0], small->p, B->Dimension);
    Matrix_Free(basis);
    Polyhedron_Free(B);

    T = unimodular_complete(small);
    Vector_Free(small);
    Vector_Exchange(T->p[0], T->p[2], T->NbColumns);
    Matrix_Print(stderr, P_VALUE_FMT, T);

    Polyhedron_Print(stderr, P_VALUE_FMT, C);
    P = Polyhedron_Image(C, T, options->MaxRays);
    Matrix_Free(T);
    Polyhedron_Free(C);

    Polyhedron_Print(stderr, P_VALUE_FMT, P);

    up.SetLength(2);
    up[0] = 1;
    up[1] = 0;

    S = barvinok_enumerate_scarf_series(P, 2, 2, options);
    S->print(std::cerr, 0, NULL);
    cerr << endl;

    S_shift = new gen_fun(S);
    S_divide = new gen_fun(S);
    S_divide->divide(up);
    S_divide->shift(up);

    value_init(c);

    do {
	S_shift->shift(up);
	hp = S->Hadamard_product(S_shift, options);
	S->add(mone, hp);
	delete hp;

	S_divide->shift(up);
	hp = S->Hadamard_product(S_divide, options);
	ok = hp->summate(&c);
	if (ok)
	    assert(value_zero_p(c));
	delete hp;
    } while(!ok);

    value_clear(c);

    frob = S->summate(1, options);
    frob->print(std::cout, 0, NULL);
    cout << endl;
    delete frob;

    delete S_shift;
    delete S_divide;
    delete S;

    Polyhedron_Free(P);
    barvinok_options_free(options);
}
