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

void dpoly::div(dpoly& d, mpq_t count, ZZ& sign)
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
    if (sign == -1)
	mpq_sub(count, count, c[len-1]);
    else
	mpq_add(count, count, c[len-1]);

    value_clear(tmp);
    mpq_clear(qtmp);
    for (int i = 0; i < len; ++i)
	mpq_clear(c[i]);
    delete [] c;
}

void dpoly_r::add_term(int i, int * powers, ZZ& coeff)
{
    if (coeff == 0)
	return;
    for (int k = 0; k < c[i].size(); ++k) {
	if (memcmp(c[i][k]->powers, powers, dim * sizeof(int)) == 0) {
	    c[i][k]->coeff += coeff;
	    return;
	}
    }
    dpoly_r_term *t = new dpoly_r_term;
    t->powers = new int[dim];
    memcpy(t->powers, powers, dim * sizeof(int));
    t->coeff = coeff;
    c[i].push_back(t);
}

dpoly_r::dpoly_r(int len, int dim)
{
    denom = 1;
    this->len = len;
    this->dim = dim;
    c = new vector< dpoly_r_term * > [len];
}

dpoly_r::dpoly_r(dpoly& num, int dim)
{
    denom = 1;
    len = num.coeff.length();
    c = new vector< dpoly_r_term * > [len];
    this->dim = dim;
    int powers[dim];
    memset(powers, 0, dim * sizeof(int));

    for (int i = 0; i < len; ++i) {
	ZZ coeff = num.coeff[i];
	add_term(i, powers, coeff);
    }
}

dpoly_r::dpoly_r(dpoly& num, dpoly& den, int pos, int dim)
{
    denom = 1;
    len = num.coeff.length();
    c = new vector< dpoly_r_term * > [len];
    this->dim = dim;
    int powers[dim];

    for (int i = 0; i < len; ++i) {
	ZZ coeff = num.coeff[i];
	memset(powers, 0, dim * sizeof(int));
	powers[pos] = 1;

	add_term(i, powers, coeff);

	for (int j = 1; j <= i; ++j) {
	    for (int k = 0; k < c[i-j].size(); ++k) {
		memcpy(powers, c[i-j][k]->powers, dim*sizeof(int));
		powers[pos]++;
		coeff = -den.coeff[j-1] * c[i-j][k]->coeff;
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
    c = new vector< dpoly_r_term * > [len];
    this->dim = dim;
    int powers[dim];
    ZZ coeff;

    for (int i = 0 ; i < len; ++i) {
	for (int k = 0; k < num->c[i].size(); ++k) {
	    memcpy(powers, num->c[i][k]->powers, dim*sizeof(int));
	    powers[pos]++;
	    add_term(i, powers, num->c[i][k]->coeff);
	}

	for (int j = 1; j <= i; ++j) {
	    for (int k = 0; k < c[i-j].size(); ++k) {
		memcpy(powers, c[i-j][k]->powers, dim*sizeof(int));
		powers[pos]++;
		coeff = -den.coeff[j-1] * c[i-j][k]->coeff;
		add_term(i, powers, coeff);
	    }
	}
    }
}

dpoly_r::~dpoly_r()
{
    for (int i = 0 ; i < len; ++i)
	for (int k = 0; k < c[i].size(); ++k) {
	    delete [] c[i][k]->powers;
	    delete c[i][k];
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
	for (int k = 0; k < c[i].size(); ++k) {
	    coeff = c[i][k]->coeff * inv_d;
	    rc->add_term(i, c[i][k]->powers, coeff);
	}

	for (int j = 1; j <= i; ++j) {
	    for (int k = 0; k < rc->c[i-j].size(); ++k) {
		coeff = - d.coeff[j] * rc->c[i-j][k]->coeff / d.coeff[0];
		rc->add_term(i, rc->c[i-j][k]->powers, coeff);
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
	for (int j = 0; j < c[i].size(); ++j) {
	    for (int k = 0; k < dim; ++k) {
		cerr << c[i][j]->powers[k] << " ";
	    }
	    cerr << ": " << c[i][j]->coeff << "/" << denom << endl;
	}
	cerr << endl;
    }
}
