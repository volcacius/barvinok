#include <iostream>
#include <vector>
#include "dpoly.h"

using std::cerr;
using std::endl;
using std::vector;

dpoly::dpoly(int d, ZZ& degree, int offset)
{
    coeff.SetLength(d+1);

    int min = d + offset;
    if (degree >= 0 && degree < ZZ(INIT_VAL, min))
	min = to_int(degree);

    ZZ c = ZZ(INIT_VAL, 1);
    if (!offset)
	coeff[0] = c;
    for (int i = 1; i <= min; ++i) {
	c *= (degree -i + 1);
	c /= i;
	coeff[i-offset] = c;
    }
}

void dpoly::operator *= (dpoly& f)
{
    assert(coeff.length() == f.coeff.length());
    vec_ZZ old = coeff;
    coeff = f.coeff[0] * coeff;
    for (int i = 1; i < coeff.length(); ++i)
	for (int j = 0; i+j < coeff.length(); ++j)
	    coeff[i+j] += f.coeff[i] * old[j];
}

mpq_t *dpoly::div(dpoly& d) const
{
    int len = coeff.length();
    Value tmp;
    value_init(tmp);
    mpq_t* c = new mpq_t[coeff.length()];
    mpq_t qtmp;
    mpq_init(qtmp);
    for (int i = 0; i < len; ++i) {
	mpq_init(c[i]);
	zz2value(coeff[i], tmp);
	mpq_set_z(c[i], tmp);

	for (int j = 1; j <= i; ++j) {
	    zz2value(d.coeff[j], tmp);
	    mpq_set_z(qtmp, tmp);
	    mpq_mul(qtmp, qtmp, c[i-j]);
	    mpq_sub(c[i], c[i], qtmp);
	}

	zz2value(d.coeff[0], tmp);
	mpq_set_z(qtmp, tmp);
	mpq_div(c[i], c[i], qtmp);
    }
    value_clear(tmp);
    mpq_clear(qtmp);

    return c;
}

void dpoly::clear_div(mpq_t *c) const
{
    int len = coeff.length();

    for (int i = 0; i < len; ++i)
	mpq_clear(c[i]);
    delete [] c;
}

void dpoly::div(dpoly& d, mpq_t count, ZZ& sign)
{
    int len = coeff.length();
    mpq_t *c = div(d);

    if (sign == -1)
	mpq_sub(count, count, c[len-1]);
    else
	mpq_add(count, count, c[len-1]);

    clear_div(c);
}

void dpoly::div(dpoly& d, mpq_t *count, const mpq_t& factor)
{
    int len = coeff.length();
    mpq_t *c = div(d);

    for (int i = 0; i < len; ++i) {
	mpq_mul(c[len-1 - i], c[len-1 - i], factor);
	mpq_add(count[i], count[i], c[len-1 - i]);
    }

    clear_div(c);
}

void dpoly_r::add_term(int i, const vector<int>& powers, const ZZ& coeff)
{
    if (coeff == 0)
	return;

    dpoly_r_term tmp;
    tmp.powers = powers;
    dpoly_r_term_list::iterator k = c[i].find(&tmp);
    if (k != c[i].end()) {
	(*k)->coeff += coeff;
	return;
    }
    dpoly_r_term *t = new dpoly_r_term;
    t->powers = powers;
    t->coeff = coeff;
    c[i].insert(t);
}

dpoly_r::dpoly_r(int len, int dim)
{
    denom = 1;
    this->len = len;
    this->dim = dim;
    c = new dpoly_r_term_list[len];
}

dpoly_r::dpoly_r(dpoly& num, int dim)
{
    denom = 1;
    len = num.coeff.length();
    c = new dpoly_r_term_list[len];
    this->dim = dim;
    vector<int> powers(dim, 0);

    for (int i = 0; i < len; ++i) {
	ZZ coeff = num.coeff[i];
	add_term(i, powers, coeff);
    }
}

dpoly_r::dpoly_r(dpoly& num, dpoly& den, int pos, int dim)
{
    denom = 1;
    len = num.coeff.length();
    c = new dpoly_r_term_list[len];
    this->dim = dim;
    int powers[dim];

    for (int i = 0; i < len; ++i) {
	ZZ coeff = num.coeff[i];
	vector<int> powers(dim, 0);
	powers[pos] = 1;

	add_term(i, powers, coeff);

	for (int j = 1; j <= i; ++j) {
	    dpoly_r_term_list::iterator k;
	    for (k = c[i-j].begin(); k != c[i-j].end(); ++k) {
		powers = (*k)->powers;
		powers[pos]++;
		coeff = -den.coeff[j-1] * (*k)->coeff;
		add_term(i, powers, coeff);
	    }
	}
    }
    //dump();
}

dpoly_r::dpoly_r(dpoly_r* num, dpoly& den, int pos, int dim)
{
    denom = num->denom;
    len = num->len;
    c = new dpoly_r_term_list[len];
    this->dim = dim;
    ZZ coeff;

    for (int i = 0 ; i < len; ++i) {
	dpoly_r_term_list::iterator k;
	for (k = num->c[i].begin(); k != num->c[i].end(); ++k) {
	    vector<int> powers = (*k)->powers;
	    powers[pos]++;
	    add_term(i, powers, (*k)->coeff);
	}

	for (int j = 1; j <= i; ++j) {
	    dpoly_r_term_list::iterator k;
	    for (k = c[i-j].begin(); k != c[i-j].end(); ++k) {
		vector<int> powers = (*k)->powers;
		powers[pos]++;
		coeff = -den.coeff[j-1] * (*k)->coeff;
		add_term(i, powers, coeff);
	    }
	}
    }
}

dpoly_r::~dpoly_r()
{
    for (int i = 0 ; i < len; ++i)
	for (dpoly_r_term_list::iterator k = c[i].begin(); k != c[i].end(); ++k) {
	    delete (*k);
	}
    delete [] c;
}

dpoly_r *dpoly_r::div(dpoly& d)
{
    dpoly_r *rc = new dpoly_r(len, dim);
    rc->denom = power(d.coeff[0], len);
    ZZ inv_d = rc->denom / d.coeff[0];
    ZZ coeff;

    for (int i = 0; i < len; ++i) {
	for (dpoly_r_term_list::iterator k = c[i].begin(); k != c[i].end(); ++k) {
	    coeff = (*k)->coeff * inv_d;
	    rc->add_term(i, (*k)->powers, coeff);
	}

	for (int j = 1; j <= i; ++j) {
	    dpoly_r_term_list::iterator k;
	    for (k = rc->c[i-j].begin(); k != rc->c[i-j].end(); ++k) {
		coeff = - d.coeff[j] * (*k)->coeff / d.coeff[0];
		rc->add_term(i, (*k)->powers, coeff);
	    }
	}
    }
    return rc;
}

void dpoly_r::dump(void)
{
    for (int i = 0; i < len; ++i) {
	cerr << endl;
	cerr << i << endl;
	cerr << c[i].size() << endl;
	for (dpoly_r_term_list::iterator j = c[i].begin(); j != c[i].end(); ++j) {
	    for (int k = 0; k < dim; ++k) {
		cerr << (*j)->powers[k] << " ";
	    }
	    cerr << ": " << (*j)->coeff << "/" << denom << endl;
	}
	cerr << endl;
    }
}
