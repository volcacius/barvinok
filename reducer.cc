#include <vector>
#include <barvinok/util.h>
#include "reducer.h"
#include "lattice_point.h"

using std::vector;

struct OrthogonalException Orthogonal;

void np_base::handle(const signed_cone& sc, barvinok_options *options)
{
    assert(sc.rays.NumRows() == dim);
    factor.n *= sc.sign;
    handle(sc.rays, current_vertex, factor, sc.closed, options);
    factor.n *= sc.sign;
}

void np_base::start(Polyhedron *P, barvinok_options *options)
{
    QQ factor(1, 1);
    for (;;) {
	try {
	    init(P);
	    for (int i = 0; i < P->NbRays; ++i) {
		if (!value_pos_p(P->Ray[i][dim+1]))
		    continue;

		Polyhedron *C = supporting_cone(P, i);
		do_vertex_cone(factor, C, P->Ray[i]+1, options);
	    }
	    break;
	} catch (OrthogonalException &e) {
	    reset();
	}
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
void normalize(ZZ& sign, vec_ZZ& num_s, mat_ZZ& num_p, vec_ZZ& den_s, vec_ZZ& den_p,
	       mat_ZZ& f)
{
    unsigned dim = f.NumRows();
    unsigned nparam = num_p.NumCols();
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
		for (int i = 0; i < num_p.NumRows(); ++i)
		    num_p[i] += f[j];
	    }
	} else
	    den_p[j] = 0;
	if (den_s[j] > 0)
	    change ^= 1;
	else {
	    den_s[j] = abs(den_s[j]);
	    for (int i = 0; i < num_p.NumRows(); ++i)
		num_s[i] += den_s[j];
	}
    }

    if (change)
	sign = -sign;
}

void reducer::base(const vec_QQ& c, const mat_ZZ& num, const mat_ZZ& den_f)
{
    for (int i = 0; i < num.NumRows(); ++i)
	base(c[i], num[i], den_f);
}

void reducer::reduce(const vec_QQ& c, const mat_ZZ& num, const mat_ZZ& den_f)
{
    assert(c.length() == num.NumRows());
    unsigned len = den_f.NumRows();  // number of factors in den

    if (num.NumCols() == lower) {
	base(c, num, den_f);
	return;
    }
    assert(num.NumCols() > 1);
    assert(num.NumRows() > 0);

    vec_ZZ den_s;
    mat_ZZ den_r;
    vec_ZZ num_s;
    mat_ZZ num_p;

    split(num, num_s, num_p, den_f, den_s, den_r);

    vec_ZZ den_p;
    den_p.SetLength(len);

    ZZ sign(INIT_VAL, 1);
    normalize(sign, num_s, num_p, den_s, den_p, den_r);

    int only_param = 0;	    // k-r-s from text
    int no_param = 0;	    // r from text
    for (int k = 0; k < len; ++k) {
	if (den_p[k] == 0)
	    ++no_param;
	else if (den_s[k] == 0)
	    ++only_param;
    }
    if (no_param == 0) {
	vec_QQ c2 = c;
	c2 *= sign;
	reduce(c2, num_p, den_r);
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

	dpoly *n[num_s.length()];
	/* We can further optimize the computation by combining the
	 * n's for rows with the same num_p
	 */
	for (int i = 0; i < num_s.length(); ++i)
	    n[i] = new dpoly(no_param, num_s[i]);
	dpoly D(no_param, den_s[k], 1);
	for ( ; ++k < len; )
	    if (den_p[k] == 0) {
		dpoly fact(no_param, den_s[k], 1);
		D *= fact;
	    }

	if (no_param + only_param == len) {
	    vec_QQ q;
	    q.SetLength(num_s.length());
	    for (int i = 0; i < num_s.length(); ++i) {
		mpq_set_si(tcount, 0, 1);
		n[i]->div(D, tcount, one);

		value2zz(mpq_numref(tcount), q[i].n);
		value2zz(mpq_denref(tcount), q[i].d);
		q[i] *= c[i];
		q[i] *= sign;
	    }
	    for (int i = q.length()-1; i >= 0; --i) {
		if (q[i].n == 0) {
		    q[i] = q[q.length()-1];
		    num_p[i] = num_p[q.length()-1];
		    q.SetLength(q.length()-1);
		    num_p.SetDims(num_p.NumRows()-1, num_p.NumCols());
		}
	    }

	    if (q.length() != 0)
		reduce(q, num_p, pden);
	} else {
	    dpoly_r *r[num_s.length()];
	    for (int i = 0; i < num_s.length(); ++i)
		r[i] = NULL;

	    for (k = 0; k < len; ++k) {
		if (den_s[k] == 0 || den_p[k] == 0)
		    continue;

		dpoly pd(no_param-1, den_s[k], 1);

		int l;
		for (l = 0; l < k; ++l)
		    if (den_r[l] == den_r[k])
			break;

		if (!r[0]) {
		    for (int i = 0; i < num_s.length(); ++i)
			r[i] = new dpoly_r(*n[i], pd, l, len);
		} else {
		    for (int i = 0; i < num_s.length(); ++i) {
			dpoly_r *nr = new dpoly_r(r[i], pd, l, len);
			delete r[i];
			r[i] = nr;
		    }
		}
	    }

	    /* We need to further optimize this part to exploit the common
	     * denominator.
	     */
	    vec_QQ factor;
	    factor.SetLength(1);
	    mat_ZZ num_p_i;
	    num_p_i.SetDims(1, num_p.NumCols());
	    int common = pden.NumRows();
	    for (int i = 0; i < num_s.length(); ++i) {
		dpoly_r *rc = r[i]->div(D);
		num_p_i[0] = num_p[i];

		factor[0].d = rc->denom * c[i].d;

		dpoly_r_term_list& final = rc->c[rc->len-1];
		int rows;
		dpoly_r_term_list::iterator j;
		for (j = final.begin(); j != final.end(); ++j) {
		    if ((*j)->coeff == 0)
			continue;
		    rows = common;
		    pden.SetDims(rows, pden.NumCols());
		    for (int k = 0; k < rc->dim; ++k) {
			int n = (*j)->powers[k];
			if (n == 0)
			    continue;
			pden.SetDims(rows+n, pden.NumCols());
			for (int l = 0; l < n; ++l)
			    pden[rows+l] = den_r[k];
			rows += n;
		    }
		    factor[0].n = (*j)->coeff *= c[i].n * sign;
		    reduce(factor, num_p_i, pden);
		}

		delete rc;
	    }

	    for (int i = 0; i < num_s.length(); ++i)
		delete r[i];
	}
	for (int i = 0; i < num_s.length(); ++i)
	    delete n[i];
    }
}

void reducer::handle(const mat_ZZ& den, Value *V, QQ c, int *closed,
		     barvinok_options *options)
{
    vec_QQ vc;
    vc.SetLength(1);
    vc[0] = c;

    lattice_point(V, den, vertex[0], closed);

    reduce(vc, vertex, den);
}

void split_one(const mat_ZZ& num, vec_ZZ& num_s, mat_ZZ& num_p,
	       const mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r)
{
    unsigned len = den_f.NumRows();  // number of factors in den
    unsigned d = num.NumCols() - 1;

    den_s.SetLength(len);
    den_r.SetDims(len, d);

    for (int r = 0; r < len; ++r) {
	den_s[r] = den_f[r][0];
	for (int k = 1; k <= d; ++k)
	    den_r[r][k-1] = den_f[r][k];
    }

    num_s.SetLength(num.NumRows());
    num_p.SetDims(num.NumRows(), d);
    for (int i = 0; i < num.NumRows(); ++i) {
	num_s[i] = num[i][0];
	for (int k = 1 ; k <= d; ++k)
	    num_p[i][k-1] = num[i][k];
    }
}

void normalize(ZZ& sign, ZZ& num, vec_ZZ& den)
{
    unsigned dim = den.length();

    int change = 0;

    for (int j = 0; j < den.length(); ++j) {
	if (den[j] > 0)
	    change ^= 1;
	else {
	    den[j] = abs(den[j]);
	    num += den[j];
	}
    }
    if (change)
	sign = -sign;
}

void icounter::base(const QQ& c, const vec_ZZ& num, const mat_ZZ& den_f)
{
    int r;
    unsigned len = den_f.NumRows();  // number of factors in den
    vec_ZZ den_s;
    den_s.SetLength(len);
    ZZ num_s = num[0];
    for (r = 0; r < len; ++r)
	den_s[r] = den_f[r][0];
    ZZ sign = ZZ(INIT_VAL, 1);
    normalize(sign, num_s, den_s);

    dpoly n(len, num_s);
    dpoly D(len, den_s[0], 1);
    for (int k = 1; k < len; ++k) {
	dpoly fact(len, den_s[k], 1);
	D *= fact;
    }
    mpq_set_si(tcount, 0, 1);
    n.div(D, tcount, one);
    zz2value(c.n, tn);
    if (sign == -1)
	value_oppose(tn, tn);
    zz2value(c.d, td);
    mpz_mul(mpq_numref(tcount), mpq_numref(tcount), tn);
    mpz_mul(mpq_denref(tcount), mpq_denref(tcount), td);
    mpq_canonicalize(tcount);
    mpq_add(count, count, tcount);
}

void infinite_icounter::base(const QQ& c, const vec_ZZ& num, const mat_ZZ& den_f)
{
    int r;
    unsigned len = den_f.NumRows();  // number of factors in den
    vec_ZZ den_s;
    den_s.SetLength(len);
    ZZ num_s = num[0];

    for (r = 0; r < len; ++r)
	den_s[r] = den_f[r][0];
    ZZ sign = ZZ(INIT_VAL, 1);
    normalize(sign, num_s, den_s);

    dpoly n(len, num_s);
    dpoly D(len, den_s[0], 1);
    for (int k = 1; k < len; ++k) {
	dpoly fact(len, den_s[k], 1);
	D *= fact;
    }

    Value tmp;
    mpq_t factor;
    mpq_init(factor);
    value_init(tmp);
    zz2value(c.n, tmp);
    if (sign == -1)
	value_oppose(tmp, tmp);
    value_assign(mpq_numref(factor), tmp);
    zz2value(c.d, tmp);
    value_assign(mpq_denref(factor), tmp);

    n.div(D, count, factor);

    value_clear(tmp);
    mpq_clear(factor);
}
