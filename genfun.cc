#include <iostream>
#include <genfun.h>

using std::cout;

void gen_fun::add(ZZ& cn, ZZ& cd, vec_ZZ& num, mat_ZZ& den)
{
    short_rat * r = new short_rat;
    r->n.coeff.SetDims(1, 2);
    r->n.coeff[0][0] = cn;
    r->n.coeff[0][1] = cd;
    r->n.power.SetDims(1, num.length());
    r->n.power[0] = num;
    r->d.power = den;

    term.push_back(r);
}

static void print_power(vec_ZZ& c, vec_ZZ& p)
{
    bool first = true;

    for (int i = 0; i < p.length(); ++i) {
	if (p[i] == 0)
	    continue;
	if (first) {
	    if (c[0] == -1 && c[1] == 1)
		cout << "-";
	    else if (c[0] != 1 || c[1] != 1) {
		cout << c;
		if (c[1] != 1)
		    cout << " / " << c[1];
		cout << " * ";
	    }
	    first = false;
	} else
	    cout << " * ";
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

void gen_fun::print(void)
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
	    if (j != 0)
		cout << " + ";
	    print_power(term[i]->n.coeff[j], term[i]->n.power[j]);
	}
	cout << ")/(";
	for (int j = 0; j < term[i]->d.power.NumRows(); ++j) {
	    if (j != 0)
		cout << " * ";
	    cout << "(1";
	    print_power(mone, term[i]->d.power[j]);
	    cout << ")";
	}
	cout << ")";
    }
}
