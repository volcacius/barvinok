#include <iostream>
#include <assert.h>
#include <genfun.h>
#include <barvinok.h>
#include <barvinok2.h>

using std::cout;

static int lex_cmp(vec_ZZ& a, vec_ZZ& b)
{
    assert(a.length() == b.length());

    for (int j = 0; j < a.length(); ++j)
	if (a[j] != b[j])
	    return a[j] < b[j] ? -1 : 1;
    return 0;
}

static int lex_cmp(mat_ZZ& a, mat_ZZ& b)
{
    assert(a.NumCols() == b.NumCols());
    int alen = a.NumRows();
    int blen = b.NumRows();
    int len = alen < blen ? alen : blen;

    for (int i = 0; i < len; ++i) {
	int s = lex_cmp(a[i], b[i]);
	if (s)
	    return s;
    }
    return alen-blen;
}

void gen_fun::add(ZZ& cn, ZZ& cd, vec_ZZ& num, mat_ZZ& den)
{
    if (cn == 0)
	return;

    short_rat * r = new short_rat;
    r->n.coeff.SetDims(1, 2);
    r->n.coeff[0][0] = cn;
    r->n.coeff[0][1] = cd;
    r->n.power.SetDims(1, num.length());
    r->n.power[0] = num;
    r->d.power = den;

    for (int i = 0; i < r->d.power.NumRows(); ++i) {
	int j;
	for (j = 0; j < r->d.power.NumCols(); ++j)
	    if (r->d.power[i][j] != 0)
		break;
	if (r->d.power[i][j] < 0) {
	    r->d.power[i] = -r->d.power[i];
	    r->n.coeff[0][0] = -r->n.coeff[0][0];
	    r->n.power[0] += r->d.power[i];
	}
    }

    for (int i = 0; i < term.size(); ++i)
	if (lex_cmp(term[i]->d.power, r->d.power) == 0) {
	    int len = term[i]->n.coeff.NumRows();
	    int dim = term[i]->n.power.NumCols();
	    term[i]->n.coeff.SetDims(len+1, 2);
	    term[i]->n.power.SetDims(len+1, dim);
	    term[i]->n.coeff[len] = r->n.coeff[0];
	    term[i]->n.power[len] = r->n.power[0];
	    delete r;
	    return;
	}

    term.push_back(r);
}

static void print_power(vec_ZZ& c, vec_ZZ& p,
			unsigned int nparam, char **param_name)
{
    bool first = true;

    for (int i = 0; i < p.length(); ++i) {
	if (p[i] == 0)
	    continue;
	if (first) {
	    if (c[0] == -1 && c[1] == 1)
		cout << "-";
	    else if (c[0] != 1 || c[1] != 1) {
		cout << c[0];
		if (c[1] != 1)
		    cout << " / " << c[1];
		cout << " * ";
	    }
	    first = false;
	} else
	    cout << " * ";
	if (i < nparam)
	    cout << param_name[i];
	else
	    cout << "x" << i;
	if (p[i] == 1)
	    continue;
	if (p[i] < 0)
	    cout << "^(" << p[i] << ")";
	else
	    cout << "^" << p[i];
    }
    if (first) {
	cout << c[0];
	if (c[1] != 1)
	    cout << " / " << c[1];
    }
}

void gen_fun::print(unsigned int nparam, char **param_name)
{
    vec_ZZ mone;
    mone.SetLength(1);
    mone[0] = -1;
    mone[1] = 1;
    for (int i = 0; i < term.size(); ++i) {
	if (i != 0)
	    cout << " + ";
	cout << "(";
	for (int j = 0; j < term[i]->n.coeff.NumRows(); ++j) {
	    if (j != 0 && term[i]->n.coeff[j][0] > 0)
		cout << "+";
	    print_power(term[i]->n.coeff[j], term[i]->n.power[j],
			nparam, param_name);
	}
	cout << ")/(";
	for (int j = 0; j < term[i]->d.power.NumRows(); ++j) {
	    if (j != 0)
		cout << " * ";
	    cout << "(1";
	    print_power(mone, term[i]->d.power[j],
			nparam, param_name);
	    cout << ")";
	}
	cout << ")";
    }
}

gen_fun::operator evalue *()
{
    evalue *EP = NULL;
    evalue factor;
    value_init(factor.d);
    value_init(factor.x.n);
    for (int i = 0; i < term.size(); ++i) {
	unsigned nvar = term[i]->d.power.NumRows();
	unsigned nparam = term[i]->d.power.NumCols();
	Matrix *C = Matrix_Alloc(nparam + nvar, 1 + nvar + nparam + 1); 
	mat_ZZ& d = term[i]->d.power;
	Polyhedron *U = Universe_Polyhedron(nparam);

	for (int j = 0; j < term[i]->n.coeff.NumRows(); ++j) {
	    for (int r = 0; r < nparam; ++r) {
		value_set_si(C->p[r][0], 0);
		for (int c = 0; c < nvar; ++c) {
		    zz2value(d[c][r], C->p[r][1+c]);
		}
		Vector_Set(&C->p[r][1+nvar], 0, nparam);
		value_set_si(C->p[r][1+nvar+r], -1);
		zz2value(term[i]->n.power[j][r], C->p[r][1+nvar+nparam]);
	    }
	    for (int r = 0; r < nvar; ++r) {
		value_set_si(C->p[nparam+r][0], 1);
		Vector_Set(&C->p[nparam+r][1], 0, nvar + nparam + 1);
		value_set_si(C->p[nparam+r][1+r], 1);
	    }
	    Polyhedron *P = Constraints2Polyhedron(C, 0);
	    evalue *E = barvinok_enumerate_ev(P, U, 0);
	    Polyhedron_Free(P);
	    zz2value(term[i]->n.coeff[j][0], factor.x.n);
	    zz2value(term[i]->n.coeff[j][1], factor.d);
	    emul(&factor, E);
	    /*
	    Matrix_Print(stdout, P_VALUE_FMT, C);
	    char *test[] = { "A", "B", "C", "D", "E", "F", "G" };
	    print_evalue(stdout, E, test);
	    */
	    if (!EP)
		EP = E;
	    else {
		eadd(E, EP);
		free_evalue_refs(E);
	    }
	}
	Matrix_Free(C);
	Polyhedron_Free(U);
    }
    value_clear(factor.d);
    value_clear(factor.x.n);
    return EP;
}
