#include "counter.h"
#include "lattice_point.h"

void counter::add_falling_powers(dpoly& n, Value degree)
{
    value_increment(n.coeff->p[0], n.coeff->p[0]);
    if (n.coeff->Size == 1)
	return;

    int min = n.coeff->Size-1;
    if (value_posz_p(degree) && value_cmp_si(degree, min) < 0)
	min = VALUE_TO_INT(degree);

    Value tmp;
    value_init(tmp);
    value_assign(tmp, degree);
    value_addto(n.coeff->p[1], n.coeff->p[1], tmp);
    for (int i = 2; i <= min; ++i) {
	value_decrement(degree, degree);
	value_multiply(tmp, tmp, degree);
	mpz_divexact_ui(tmp, tmp, i);
	value_addto(n.coeff->p[i], n.coeff->p[i], tmp);
    }
    value_clear(tmp);
}

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

    dpoly d(dim);
    for (int k = 0; k < num->Size; ++k)
	add_falling_powers(d, num->p[k]);
    dpoly n(dim, den->p[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den->p[k], 1);
	n *= fact;
    }
    d.div(n, count, sign);
}
