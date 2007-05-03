#include "counter.h"
#include "lattice_point.h"

void counter::handle(const mat_ZZ& rays, Value *V, const QQ& c, unsigned long det,
		     int *closed, barvinok_options *options)
{
    for (int k = 0; k < dim; ++k) {
	if (lambda * rays[k] == 0)
	    throw Orthogonal;
    }

    assert(c.d == 1);
    assert(c.n == 1 || c.n == -1);
    sign = c.n;

    lattice_point(V, rays, vertex, det, closed);
    num = vertex * lambda;
    den = rays * lambda;
    offset = 0;
    normalize(sign, offset, den);

    num[0] += offset;
    zz2value(num[0], tz);
    dpoly d(dim, tz);
    for (int k = 1; k < num.length(); ++k) {
	num[k] += offset;
	zz2value(num[k], tz);
	dpoly term(dim, tz);
	d += term;
    }
    zz2value(den[0], tz);
    dpoly n(dim, tz, 1);
    for (int k = 1; k < dim; ++k) {
	zz2value(den[k], tz);
	dpoly fact(dim, tz, 1);
	n *= fact;
    }
    d.div(n, count, sign);
}
