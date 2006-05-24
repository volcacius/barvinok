#include <iostream>
#include <assert.h>
#include "config.h"
#include <barvinok/genfun.h>
#include <barvinok/barvinok.h>

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

void gen_fun::add(const ZZ& cn, const ZZ& cd, const vec_ZZ& num, 
		  const mat_ZZ& den)
{
    if (cn == 0)
	return;

    short_rat * r = new short_rat;
    r->n.coeff.SetDims(1, 2);
    ZZ g = GCD(cn, cd);
    r->n.coeff[0][0] = cn/g;
    r->n.coeff[0][1] = cd/g;
    r->n.power.SetDims(1, num.length());
    r->n.power[0] = num;
    r->d.power = den;

    /* Make all powers in denominator lexico-positive */
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

    /* Order powers in denominator */
    for (int i = 0; i < r->d.power.NumRows(); ++i) {
	int m = i;
	for (int j = i+1; j < r->d.power.NumRows(); ++j)
	    if (lex_cmp(r->d.power[j], r->d.power[m]) < 0)
		m = j;
	if (m != i) {
	    vec_ZZ tmp = r->d.power[m];
	    r->d.power[m] = r->d.power[i];
	    r->d.power[i] = tmp;
	}
    }

    for (int i = 0; i < term.size(); ++i)
	if (lex_cmp(term[i]->d.power, r->d.power) == 0) {
	    int len = term[i]->n.coeff.NumRows();
	    int j;
	    for (j = 0; j < len; ++j)
		if (r->n.power[0] == term[i]->n.power[j])
		    break;
	    if (j < len) {
		ZZ g = GCD(r->n.coeff[0][1], term[i]->n.coeff[j][1]);
		ZZ n = term[i]->n.coeff[j][0] * (r->n.coeff[0][1] / g) +
			   (term[i]->n.coeff[j][1] / g) * r->n.coeff[0][0];
		ZZ d = term[i]->n.coeff[j][1] / g * r->n.coeff[0][1];
		if (n != 0) {
		    g = GCD(n,d);
		    term[i]->n.coeff[j][0] = n/g;
		    term[i]->n.coeff[j][1] = d/g;
		} else {
		    if (len > 1) {
			if (j < len-1) {
			    term[i]->n.power[j] = term[i]->n.power[len-1];
			    term[i]->n.coeff[j] = term[i]->n.coeff[len-1];
			}
			int dim = term[i]->n.power.NumCols();
			term[i]->n.coeff.SetDims(len-1, 2);
			term[i]->n.power.SetDims(len-1, dim);
		    } else {
			delete term[i];
			if (i != term.size()-1)
			    term[i] = term[term.size()-1];
			term.pop_back();
		    }
		}
	    } else {
		int dim = term[i]->n.power.NumCols();
		term[i]->n.coeff.SetDims(len+1, 2);
		term[i]->n.power.SetDims(len+1, dim);
		term[i]->n.coeff[len] = r->n.coeff[0];
		term[i]->n.power[len] = r->n.power[0];
	    }
	    delete r;
	    return;
	}

    term.push_back(r);
}

void gen_fun::add(const ZZ& cn, const ZZ& cd, gen_fun *gf)
{
    ZZ n, d;
    for (int i = 0; i < gf->term.size(); ++i) {
	for (int j = 0; j < gf->term[i]->n.power.NumRows(); ++j) {
	    n = cn * gf->term[i]->n.coeff[j][0];
	    d = cd * gf->term[i]->n.coeff[j][1];
	    add(n, d, gf->term[i]->n.power[j], gf->term[i]->d.power);
	}
    }
}

/*
 * Perform the substitution specified by CP and (map, offset)
 *
 * CP is a homogeneous matrix that maps a set of "compressed parameters"
 * to the original set of parameters.
 *
 * This function is applied to a gen_fun computed with the compressed parameters
 * and adapts it to refer to the original parameters.
 *
 * That is, if y are the compressed parameters and x = A y + b are the original
 * parameters, then we want the coefficient of the monomial t^y in the original
 * generating function to be the coefficient of the monomial u^x in the resulting
 * generating function.  
 * The original generating function has the form
 *
 *        a t^m/(1-t^n) = a t^m + a t^{m+n} + a t^{m+2n} + ...
 * 
 * Since each term t^y should correspond to a term u^x, with x = A y + b, we want
 *
 *         a u^{A m + b} + a u^{A (m+n) + b} + a u^{A (m+2n) +b} + ... = 
 *        
 *         = a u^{A m + b}/(1-u^{A n})
 *
 * Therefore, we multiply the powers m and n in both numerator and denominator by A
 * and add b to the power in the numerator.
 * Since the above powers are stored as row vectors m^T and n^T,
 * we compute, say, m'^T = m^T A^T to obtain m' = A m.
 *
 * The pair (map, offset) contains the same information as CP.
 * map is the transpose of the linear part of CP, while offset is the constant part.
 */
void gen_fun::substitute(Matrix *CP, const mat_ZZ& map, const vec_ZZ& offset)
{
    Polyhedron *C = Polyhedron_Image(context, CP, 0);
    Polyhedron_Free(context);
    context = C;
    for (int i = 0; i < term.size(); ++i) {
	term[i]->d.power *= map;
	term[i]->n.power *= map;
	for (int j = 0; j < term[i]->n.power.NumRows(); ++j)
	    term[i]->n.power[j] += offset;
    }
}

gen_fun *gen_fun::Hadamard_product(gen_fun *gf, unsigned MaxRays)
{
    Polyhedron *C = DomainIntersection(context, gf->context, MaxRays);
    Polyhedron *U = Universe_Polyhedron(C->Dimension);
    gen_fun *sum = new gen_fun(C);
    for (int i = 0; i < term.size(); ++i) {
	for (int i2 = 0; i2 < gf->term.size(); ++i2) {
	    int d = term[i]->d.power.NumCols();
	    int k1 = term[i]->d.power.NumRows();
	    int k2 = gf->term[i2]->d.power.NumRows();
	    assert(term[i]->d.power.NumCols() == gf->term[i2]->d.power.NumCols());
	    for (int j = 0; j < term[i]->n.power.NumRows(); ++j) {
		for (int j2 = 0; j2 < gf->term[i2]->n.power.NumRows(); ++j2) {
		    Matrix *M = Matrix_Alloc(k1+k2+d+d, 1+k1+k2+d+1);
		    for (int k = 0; k < k1+k2; ++k) {
			value_set_si(M->p[k][0], 1);
			value_set_si(M->p[k][1+k], 1);
		    }
		    for (int k = 0; k < d; ++k) {
			value_set_si(M->p[k1+k2+k][1+k1+k2+k], -1);
			zz2value(term[i]->n.power[j][k], M->p[k1+k2+k][1+k1+k2+d]);
			for (int l = 0; l < k1; ++l)
			    zz2value(term[i]->d.power[l][k], M->p[k1+k2+k][1+l]);
		    }
		    for (int k = 0; k < d; ++k) {
			value_set_si(M->p[k1+k2+d+k][1+k1+k2+k], -1);
			zz2value(gf->term[i2]->n.power[j2][k], 
				 M->p[k1+k2+d+k][1+k1+k2+d]);
			for (int l = 0; l < k2; ++l)
			    zz2value(gf->term[i2]->d.power[l][k], 
				     M->p[k1+k2+d+k][1+k1+l]);
		    }
		    Polyhedron *P = Constraints2Polyhedron(M, MaxRays);
		    Matrix_Free(M);

		    gen_fun *t = barvinok_series(P, U, MaxRays);

		    ZZ cn = term[i]->n.coeff[j][0] * gf->term[i2]->n.coeff[j2][0];
		    ZZ cd = term[i]->n.coeff[j][1] * gf->term[i2]->n.coeff[j2][1];
		    sum->add(cn, cd, t);
		    delete t;

		    Polyhedron_Free(P);
		}
	    }
	}
    }
    Polyhedron_Free(U);
    return sum;
}

void gen_fun::add_union(gen_fun *gf, unsigned MaxRays)
{
    ZZ one, mone;
    one = 1;
    mone = -1;

    gen_fun *hp = gf->Hadamard_product(gf, MaxRays);
    add(one, one, gf);
    add(mone, one, hp);
    delete hp;
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
		cout << "*";
	    }
	    first = false;
	} else
	    cout << "*";
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

void gen_fun::print(unsigned int nparam, char **param_name) const
{
    vec_ZZ mone;
    mone.SetLength(2);
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

gen_fun::operator evalue *() const
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
	Polyhedron *U = context ? context : Universe_Polyhedron(nparam);

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
	    if (EVALUE_IS_ZERO(*E)) {
		free_evalue_refs(E);
		free(E);
		continue;
	    }
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
		free(E);
	    }
	}
	Matrix_Free(C);
	if (!context)
	    Polyhedron_Free(U);
    }
    value_clear(factor.d);
    value_clear(factor.x.n);
    return EP;
}

void gen_fun::coefficient(Value* params, Value* c) const
{
    if (context && !in_domain(context, params)) {
	value_set_si(*c, 0);
	return;
    }

    evalue part;
    value_init(part.d);
    value_init(part.x.n);
    evalue sum;
    value_init(sum.d);
    evalue_set_si(&sum, 0, 1);
    Value tmp;
    value_init(tmp);

    for (int i = 0; i < term.size(); ++i) {
	unsigned nvar = term[i]->d.power.NumRows();
	unsigned nparam = term[i]->d.power.NumCols();
	Matrix *C = Matrix_Alloc(nparam + nvar, 1 + nvar + 1); 
	mat_ZZ& d = term[i]->d.power;

	for (int j = 0; j < term[i]->n.coeff.NumRows(); ++j) {
	    for (int r = 0; r < nparam; ++r) {
		value_set_si(C->p[r][0], 0);
		for (int c = 0; c < nvar; ++c) {
		    zz2value(d[c][r], C->p[r][1+c]);
		}
		zz2value(term[i]->n.power[j][r], C->p[r][1+nvar]);
		value_subtract(C->p[r][1+nvar], C->p[r][1+nvar], params[r]);
	    }
	    for (int r = 0; r < nvar; ++r) {
		value_set_si(C->p[nparam+r][0], 1);
		Vector_Set(&C->p[nparam+r][1], 0, nvar + 1);
		value_set_si(C->p[nparam+r][1+r], 1);
	    }
	    Polyhedron *P = Constraints2Polyhedron(C, 0);
	    if (emptyQ(P)) {
		Polyhedron_Free(P);
		continue;
	    }
	    barvinok_count(P, &tmp, 0);
	    Polyhedron_Free(P);
	    if (value_zero_p(tmp))
		continue;
	    zz2value(term[i]->n.coeff[j][0], part.x.n);
	    zz2value(term[i]->n.coeff[j][1], part.d);
	    value_multiply(part.x.n, part.x.n, tmp);
	    eadd(&part, &sum);
	}
	Matrix_Free(C);
    }

    assert(value_one_p(sum.d));
    value_assign(*c, sum.x.n);

    value_clear(tmp);
    value_clear(part.d);
    value_clear(part.x.n);
    value_clear(sum.d);
    value_clear(sum.x.n);
}
