#include <assert.h>
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

    assert(c.d == 1);
    assert(c.n == 1 || c.n == -1);
    int sign = to_int(c.n);

    Matrix_Vector_Product(Rays, lambda->p, den->p_Init);
    for (int k = 0; k < dim; ++k)
	if (value_zero_p(den->p_Init[k])) {
	    Matrix_Free(Rays);
	    throw Orthogonal;
	}
    Inner_Product(lambda->p, V, dim, &tmp);
    lattice_points_fixed(V, &tmp, Rays, den, num, det, closed);
    num->NbRows = det;
    Matrix_Free(Rays);

    if (dim % 2)
	sign = -sign;

    dpoly d(dim);
    for (int k = 0; k < num->NbRows; ++k)
	add_falling_powers(d, num->p_Init[k]);
    dpoly n(dim, den->p_Init[0], 1);
    for (int k = 1; k < dim; ++k) {
	dpoly fact(dim, den->p_Init[k], 1);
	n *= fact;
    }
    d.div(n, count, sign);
}
