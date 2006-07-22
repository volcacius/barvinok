#include "bfcounter.h"

static int lex_cmp(vec_ZZ& a, vec_ZZ& b)
{
    assert(a.length() == b.length());

    for (int j = 0; j < a.length(); ++j)
	if (a[j] != b[j])
	    return a[j] < b[j] ? -1 : 1;
    return 0;
}

void bf_base::add_term(bfc_term_base *t, vec_ZZ& num_orig, vec_ZZ& extra_num)
{
    vec_ZZ num;
    int d = num_orig.length();
    num.SetLength(d-1);
    for (int l = 0; l < d-1; ++l)
	num[l] = num_orig[l+1] + extra_num[l];

    add_term(t, num);
}

void bf_base::add_term(bfc_term_base *t, vec_ZZ& num)
{
    int len = t->terms.NumRows();
    int i, r;
    for (i = 0; i < len; ++i) {
	r = lex_cmp(t->terms[i], num);
	if (r >= 0)
	    break;
    }
    if (i == len || r > 0) {
	t->terms.SetDims(len+1, num.length());
	insert_term(t, i);
	t->terms[i] = num;
    } else {
	// i < len && r == 0
	update_term(t, i);
    }
}

bfc_term_base* bfcounter_base::new_bf_term(int len)
{
    bfc_term* t = new bfc_term(len);
    t->c.SetLength(0);
    return t;
}

void bfcounter_base::set_factor(bfc_term_base *t, int k, int change)
{
    bfc_term* bfct = static_cast<bfc_term *>(t);
    c = bfct->c[k];
    if (change)
	c.n = -c.n;
}

void bfcounter_base::set_factor(bfc_term_base *t, int k, mpq_t &f, int change)
{
    bfc_term* bfct = static_cast<bfc_term *>(t);
    value2zz(mpq_numref(f), c.n);
    value2zz(mpq_denref(f), c.d);
    c *= bfct->c[k];
    if (change)
	c.n = -c.n;
}

void bfcounter_base::set_factor(bfc_term_base *t, int k, const QQ& c_factor,
				int change)
{
    bfc_term* bfct = static_cast<bfc_term *>(t);
    c = bfct->c[k];
    c *= c_factor;
    if (change)
	c.n = -c.n;
}

void bfcounter_base::insert_term(bfc_term_base *t, int i)
{
    bfc_term* bfct = static_cast<bfc_term *>(t);
    int len = t->terms.NumRows()-1;	// already increased by one

    bfct->c.SetLength(len+1);
    for (int j = len; j > i; --j) {
	bfct->c[j] = bfct->c[j-1];
	t->terms[j] = t->terms[j-1];
    }
    bfct->c[i] = c;
}

void bfcounter_base::update_term(bfc_term_base *t, int i)
{
    bfc_term* bfct = static_cast<bfc_term *>(t);

    bfct->c[i] += c;
}
