#include <barvinok/util.h>
#include "tcounter.h"
#include "lattice_point.h"

using std::cerr;
using std::endl;

void tcounter::setup_todd(unsigned dim)
{
    value_set_si(todd.coeff->p[0], 1);

    dpoly d(dim);
    value_set_si(d.coeff->p[dim], 1);
    for (int i = dim-1; i >= 0; --i)
	mpz_mul_ui(d.coeff->p[i], d.coeff->p[i+1], i+2);

    todd_denom = todd.div(d);
    /* shift denominators up -> divide by (dim+1)! */
    for (int i = todd.coeff->Size-1; i >= 1; --i)
	value_assign(todd_denom->p[i], todd_denom->p[i-1]);
    value_set_si(todd_denom->p[0], 1);
}

void tcounter::adapt_todd(dpoly& t, const Value c)
{
    if (t.coeff->Size <= 1)
	return;
    value_assign(tmp, c);
    value_multiply(t.coeff->p[1], t.coeff->p[1], tmp);
    for (int i = 2; i < t.coeff->Size; ++i) {
	value_multiply(tmp, tmp, c);
	value_multiply(t.coeff->p[i], t.coeff->p[i], tmp);
    }
}

void tcounter::add_powers(dpoly& n, const Value c)
{
    value_increment(n.coeff->p[0], n.coeff->p[0]);
    if (n.coeff->Size == 1)
	return;
    value_assign(tmp, c);
    value_addto(n.coeff->p[1], n.coeff->p[1], tmp);
    for (int i = 2; i < n.coeff->Size; ++i) {
	value_multiply(tmp, tmp, c);
	value_addto(n.coeff->p[i], n.coeff->p[i], tmp);
    }
}

void tcounter::handle(const mat_ZZ& rays, Value *V, const QQ& c, unsigned long det,
		     int *closed, barvinok_options *options)
{
    for (int k = 0; k < dim; ++k) {
	if (lambda * rays[k] == 0)
	    throw Orthogonal;
    }

    assert(c.d == 1);
    assert(c.n == 1 || c.n == -1);
    sign = c.n;

    if (dim % 2)
	sign = -sign;

    lattice_point(V, rays, vertex, det, closed);
    num = vertex * lambda;
    den = rays * lambda;

    dpoly t(todd);
    zz2value(den[0], tz);
    value_assign(denom, tz);
    adapt_todd(t, tz);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(todd);
	zz2value(den[k], tz);
	value_multiply(denom, denom, tz);
	adapt_todd(fact, tz);
	t *= fact;
    }

    dpoly n(dim);
    for (int k = 0; k < num.length(); ++k) {
	zz2value(num[k], tz);
	add_powers(n, tz);
    }

    for (int i = 0; i < n.coeff->Size; ++i)
	value_multiply(n.coeff->p[i], n.coeff->p[i], todd_denom->p[i]);
    value_multiply(denom, denom, todd_denom->p[todd_denom->Size-1]);

    value_set_si(tmp, 1);
    for (int i = 2; i < n.coeff->Size; ++i) {
	mpz_mul_ui(tmp, tmp, i);
	mpz_divexact(n.coeff->p[i], n.coeff->p[i], tmp);
    }

    value_multiply(tmp, t.coeff->p[0], n.coeff->p[n.coeff->Size-1]);
    for (int i = 1; i < n.coeff->Size; ++i)
	value_addmul(tmp, t.coeff->p[i], n.coeff->p[n.coeff->Size-1-i]);

    value_assign(mpq_numref(tcount), tmp);
    value_assign(mpq_denref(tcount), denom);
    mpq_canonicalize(tcount);
    if (sign == -1)
	mpq_sub(count, count, tcount);
    else
	mpq_add(count, count, tcount);
}