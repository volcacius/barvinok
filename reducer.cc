#include <vector>
#include <barvinok/util.h>
#include "reducer.h"
#include "lattice_point.h"

using std::vector;

void np_base::handle_polar(Polyhedron *C, int s)
{
    assert(C->NbRays-1 == dim);
    factor.n = s;
    handle_polar(C, current_vertex, factor);
}

void np_base::start(Polyhedron *P, unsigned MaxRays)
{
    init(P);
    for (int i = 0; i < P->NbRays; ++i) {
	if (!value_pos_p(P->Ray[i][dim+1]))
	    continue;

	current_vertex = P->Ray[i]+1;
	Polyhedron *C = supporting_cone(P, i);
	decompose(C, MaxRays);
    }
}

/* input:
 *	f: the powers in the denominator for the remaining vars
 *	      each row refers to a factor
 *      den_s: for each factor, the power of  (s+1)
 *	sign
 *	num_s: powers in the numerator corresponding to the summed vars
 *	num_p: powers in the numerator corresponding to the remaining vars
 * number of rays in cone: "dim" = "k"
 * length of each ray: "dim" = "d"
 * for now, it is assumed: k == d
 * output:
 *	den_p: for each factor
 *		0: independent of remaining vars
 *		1: power corresponds to corresponding row in f
 *
 * all inputs are subject to change
 */
void normalize(ZZ& sign, ZZ& num_s, vec_ZZ& num_p, vec_ZZ& den_s, vec_ZZ& den_p,
	       mat_ZZ& f)
{
    unsigned dim = f.NumRows();
    unsigned nparam = num_p.length();
    unsigned nvar = dim - nparam;

    int change = 0;

    for (int j = 0; j < den_s.length(); ++j) {
	if (den_s[j] == 0) {
	    den_p[j] = 1;
	    continue;
	}
	int k;
	for (k = 0; k < nparam; ++k)
	    if (f[j][k] != 0)
		break;
	if (k < nparam) {
	    den_p[j] = 1;
	    if (den_s[j] > 0) {
		f[j] = -f[j];
		num_p += f[j];
	    }
	} else
	    den_p[j] = 0;
	if (den_s[j] > 0)
	    change ^= 1;
	else {
	    den_s[j] = abs(den_s[j]);
	    num_s += den_s[j];
	}
    }

    if (change)
	sign = -sign;
}

void reducer::reduce(QQ c, vec_ZZ& num, mat_ZZ& den_f)
{
    unsigned len = den_f.NumRows();  // number of factors in den

    if (num.length() == lower) {
	base(c, num, den_f);
	return;
    }
    assert(num.length() > 1);

    vec_ZZ den_s;
    mat_ZZ den_r;
    ZZ num_s;
    vec_ZZ num_p;

    split(num, num_s, num_p, den_f, den_s, den_r);

    vec_ZZ den_p;
    den_p.SetLength(len);

    normalize(c.n, num_s, num_p, den_s, den_p, den_r);

    int only_param = 0;	    // k-r-s from text
    int no_param = 0;	    // r from text
    for (int k = 0; k < len; ++k) {
	if (den_p[k] == 0)
	    ++no_param;
	else if (den_s[k] == 0)
	    ++only_param;
    }
    if (no_param == 0) {
	reduce(c, num_p, den_r);
    } else {
	int k, l;
	mat_ZZ pden;
	pden.SetDims(only_param, den_r.NumCols());

	for (k = 0, l = 0; k < len; ++k)
	    if (den_s[k] == 0)
		pden[l++] = den_r[k];

	for (k = 0; k < len; ++k)
	    if (den_p[k] == 0)
		break;

	dpoly n(no_param, num_s);
	dpoly D(no_param, den_s[k], 1);
	for ( ; ++k < len; )
	    if (den_p[k] == 0) {
		dpoly fact(no_param, den_s[k], 1);
		D *= fact;
	    }

	if (no_param + only_param == len) {
	    mpq_set_si(tcount, 0, 1);
	    n.div(D, tcount, one);

	    QQ q;
	    value2zz(mpq_numref(tcount), q.n);
	    value2zz(mpq_denref(tcount), q.d);

	    q *= c;

	    if (q.n != 0)
		reduce(q, num_p, pden);
	} else {
	    dpoly_r * r = 0;

	    for (k = 0; k < len; ++k) {
		if (den_s[k] == 0 || den_p[k] == 0)
		    continue;

		dpoly pd(no_param-1, den_s[k], 1);

		int l;
		for (l = 0; l < k; ++l)
		    if (den_r[l] == den_r[k])
			break;

		if (r == 0)
		    r = new dpoly_r(n, pd, l, len);
		else {
		    dpoly_r *nr = new dpoly_r(r, pd, l, len);
		    delete r;
		    r = nr;
		}
	    }

	    dpoly_r *rc = r->div(D);

	    QQ factor;
	    factor.d = rc->denom * c.d;

	    int common = pden.NumRows();
	    vector< dpoly_r_term * >& final = rc->c[rc->len-1];
	    int rows;
	    for (int j = 0; j < final.size(); ++j) {
		if (final[j]->coeff == 0)
		    continue;
		rows = common;
		pden.SetDims(rows, pden.NumCols());
		for (int k = 0; k < rc->dim; ++k) {
		    int n = final[j]->powers[k];
		    if (n == 0)
			continue;
		    pden.SetDims(rows+n, pden.NumCols());
		    for (int l = 0; l < n; ++l)
			pden[rows+l] = den_r[k];
		    rows += n;
		}
		factor.n = final[j]->coeff *= c.n;
		reduce(factor, num_p, pden);
	    }

	    delete rc;
	    delete r;
	}
    }
}

void reducer::handle_polar(Polyhedron *C, Value *V, QQ c)
{
    lattice_point(V, C, vertex);

    mat_ZZ den;
    den.SetDims(dim, dim);

    int r;
    for (r = 0; r < dim; ++r)
	values2zz(C->Ray[r]+1, den[r], dim);

    reduce(c, vertex, den);
}

void ireducer::split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
		     mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r)
{
    unsigned len = den_f.NumRows();  // number of factors in den
    unsigned d = num.length() - 1;

    den_s.SetLength(len);
    den_r.SetDims(len, d);

    for (int r = 0; r < len; ++r) {
	den_s[r] = den_f[r][0];
	for (int k = 1; k <= d; ++k)
	    den_r[r][k-1] = den_f[r][k];
    }

    num_s = num[0];
    num_p.SetLength(d);
    for (int k = 1 ; k <= d; ++k)
	num_p[k-1] = num[k];
}
