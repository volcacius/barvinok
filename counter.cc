#include "counter.h"
#include "lattice_point.h"

void counter::handle(const mat_ZZ& rays, Value *V, const QQ& c, unsigned long det,
		     int *closed, barvinok_options *options)
{
    Matrix* Rays = zz2matrix(rays);
    for (int k = 0; k < dim; ++k) {
	Inner_Product(lambda->p, Rays->p[k], dim, &num->p[0]);
	if (value_zero_p(num->p[0]))
	    throw Orthogonal;
    }

    assert(c.d == 1);
    assert(c.n == 1 || c.n == -1);
    sign = c.n;

    lattice_point(V, Rays, vertex, det, closed);
    vertex->NbRows = det;
    num->Size = det;
    Matrix_Vector_Product(vertex, lambda->p, num->p);
    Matrix_Vector_Product(Rays, lambda->p, den->p);
    Matrix_Free(Rays);

    if (dim % 2)
	sign = -sign;

    dpoly d(dim, num->p[0]);
    for (int k = 1; k < num->Size; ++k) {
	dpoly term(dim, num->p[k]);
	d += term;
    }
    dpoly n(dim, den->p[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den->p[k], 1);
	n *= fact;
    }
    d.div(n, count, sign);
}
