#include <iostream>
#include <vector>
#include "dpoly.h"

using std::cerr;
using std::endl;
using std::vector;

/* Construct truncated expansion of (1+t)^(degree),
 * computing the first 1+d coefficients
 */
dpoly::dpoly(int d, const Value degree, int offset)
{
    coeff = Vector_Alloc(d+1);

    /* For small degrees, we only need to compute some coefficients */
    int min = d + offset;
    if (value_posz_p(degree) && value_cmp_si(degree, min) < 0)
	min = VALUE_TO_INT(degree);

    Value c, tmp;
    value_init(c);
    value_init(tmp);
    value_set_si(c, 1);
    if (!offset)
	value_assign(coeff->p[0], c);
    value_assign(tmp, degree);
    for (int i = 1; i <= min; ++i) {
	value_multiply(c, c, tmp);
	value_decrement(tmp, tmp);
	mpz_divexact_ui(c, c, i);
	value_assign(coeff->p[i-offset], c);
    }
    value_clear(c);
    value_clear(tmp);
}

void dpoly::operator += (const dpoly& t)
{
    assert(coeff->Size == t.coeff->Size);
    for (int i = 0; i < coeff->Size; ++i)
	value_addto(coeff->p[i], coeff->p[i], t.coeff->p[i]);
}

void dpoly::operator *= (const Value f)
{
    for (int i = 0; i < coeff->Size; ++i)
	value_multiply(coeff->p[i], coeff->p[i], f);
}

void dpoly::operator *= (const dpoly& f)
{
    assert(coeff->Size == f.coeff->Size);
    Vector *old = Vector_Alloc(coeff->Size);
    Vector_Copy(coeff->p, old->p, coeff->Size);
    Vector_Scale(coeff->p, coeff->p, f.coeff->p[0], coeff->Size);
    for (int i = 1; i < coeff->Size; ++i)
	for (int j = 0; i+j < coeff->Size; ++j)
	    value_addmul(coeff->p[i+j], f.coeff->p[i], old->p[j]);
    Vector_Free(old);
}

mpq_t *dpoly::div(dpoly& d) const
{
    int len = coeff->Size;
    mpq_t* c = new mpq_t[coeff->Size];
    mpq_t qtmp;
    mpq_init(qtmp);
    for (int i = 0; i < len; ++i) {
	mpq_init(c[i]);
	mpq_set_z(c[i], coeff->p[i]);

	for (int j = 1; j <= i; ++j) {
	    mpq_set_z(qtmp, d.coeff->p[j]);
	    mpq_mul(qtmp, qtmp, c[i-j]);
	    mpq_sub(c[i], c[i], qtmp);
	}

	mpq_set_z(qtmp, d.coeff->p[0]);
	mpq_div(c[i], c[i], qtmp);
    }
    mpq_clear(qtmp);

    return c;
}

void dpoly::clear_div(mpq_t *c) const
{
    int len = coeff->Size;

    for (int i = 0; i < len; ++i)
	mpq_clear(c[i]);
    delete [] c;
}

void dpoly::div(dpoly& d, mpq_t count, ZZ& sign)
{
    int len = coeff->Size;
    mpq_t *c = div(d);

    if (sign == -1)
	mpq_sub(count, count, c[len-1]);
    else
	mpq_add(count, count, c[len-1]);

    clear_div(c);
}

void dpoly::div(dpoly& d, mpq_t *count, const mpq_t& factor)
{
    int len = coeff->Size;
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
    len = num.coeff->Size;
    c = new dpoly_r_term_list[len];
    this->dim = dim;
    vector<int> powers(dim, 0);

    for (int i = 0; i < len; ++i) {
	ZZ coeff;
	value2zz(num.coeff->p[i], coeff);
	add_term(i, powers, coeff);
    }
}

dpoly_r::dpoly_r(dpoly& num, dpoly& den, int pos, int dim)
{
    denom = 1;
    len = num.coeff->Size;
    c = new dpoly_r_term_list[len];
    this->dim = dim;
    int powers[dim];
    ZZ coeff;

    for (int i = 0; i < len; ++i) {
	vector<int> powers(dim, 0);
	powers[pos] = 1;

	value2zz(num.coeff->p[i], coeff);
	add_term(i, powers, coeff);

	for (int j = 1; j <= i; ++j) {
	    dpoly_r_term_list::iterator k;
	    for (k = c[i-j].begin(); k != c[i-j].end(); ++k) {
		powers = (*k)->powers;
		powers[pos]++;
		value2zz(den.coeff->p[j-1], coeff);
		negate(coeff, coeff);
		coeff *= (*k)->coeff;
		add_term(i, powers, coeff);
	    }
	}
    }
    //dump();
}

dpoly_r::dpoly_r(const dpoly_r* num, dpoly& den, int pos, int dim)
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
		value2zz(den.coeff->p[j-1], coeff);
		negate(coeff, coeff);
		coeff *= (*k)->coeff;
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

dpoly_r *dpoly_r::div(const dpoly& d) const
{
    dpoly_r *rc = new dpoly_r(len, dim);
    ZZ coeff;
    ZZ coeff0;
    value2zz(d.coeff->p[0], coeff0);
    rc->denom = power(coeff0, len);
    ZZ inv_d = rc->denom / coeff0;

    for (int i = 0; i < len; ++i) {
	for (dpoly_r_term_list::iterator k = c[i].begin(); k != c[i].end(); ++k) {
	    coeff = (*k)->coeff * inv_d;
	    rc->add_term(i, (*k)->powers, coeff);
	}

	for (int j = 1; j <= i; ++j) {
	    dpoly_r_term_list::iterator k;
	    for (k = rc->c[i-j].begin(); k != rc->c[i-j].end(); ++k) {
		value2zz(d.coeff->p[j], coeff);
		coeff = - coeff * (*k)->coeff / coeff0;
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
