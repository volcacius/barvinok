#include <vector>
#include <barvinok/util.h>
#include "reducer.h"
#include "lattice_point.h"

using std::vector;
using std::cerr;
using std::endl;

struct OrthogonalException Orthogonal;

void np_base::handle(const signed_cone& sc, barvinok_options *options)
{
    assert(sc.rays.NumRows() == dim);
    factor.n *= sc.sign;
    handle(sc.rays, current_vertex, factor, sc.det, options);
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

struct dpoly_r_scanner {
    const dpoly_r *rc;
    const dpoly * const *num;
    int n;
    int dim;
    dpoly_r_term_list::iterator *iter;
    vector<int> powers;
    vec_ZZ coeff;

    dpoly_r_scanner(const dpoly * const *num, int n, const dpoly_r *rc, int dim)
		    : num(num), rc(rc), n(n), dim(dim), powers(dim, 0) {
	coeff.SetLength(n);
	iter = new dpoly_r_term_list::iterator[rc->len];
	for (int i = 0; i < rc->len; ++i) {
	    int k;
	    for (k = 0; k < n; ++k)
		if (value_notzero_p(num[k]->coeff->p[rc->len-1-i]))
		    break;
	    if (k < n)
		iter[i] = rc->c[i].begin();
	    else
		iter[i] = rc->c[i].end();
	}
    }
    bool next() {
	int pos[rc->len];
	int len = 0;

	for (int i = 0; i < rc->len; ++i) {
	    if (iter[i] == rc->c[i].end())
		continue;
	    if (!len)
		pos[len++] = i;
	    else {
		if ((*iter[i])->powers < (*iter[pos[0]])->powers) {
		    pos[0] = i;
		    len = 1;
		} else if ((*iter[i])->powers == (*iter[pos[0]])->powers)
		    pos[len++] = i;
	    }
	}

	if (!len)
	    return false;

	powers = (*iter[pos[0]])->powers;
	for (int k = 0; k < n; ++k) {
	    value2zz(num[k]->coeff->p[rc->len-1-pos[0]], tmp);
	    mul(coeff[k], (*iter[pos[0]])->coeff, tmp);
	}
	++iter[pos[0]];
	for (int i = 1; i < len; ++i) {
	    for (int k = 0; k < n; ++k) {
		value2zz(num[k]->coeff->p[rc->len-1-pos[i]], tmp);
		mul(tmp, (*iter[pos[i]])->coeff, tmp);
		add(coeff[k], coeff[k], tmp);
	    }
	    ++iter[pos[i]];
	}

	return true;
    }
    ~dpoly_r_scanner() {
	delete [] iter;
    }
private:
    ZZ tmp;
};

void reducer::reduce_canonical(const vec_QQ& c, const mat_ZZ& num,
				const mat_ZZ& den_f)
{
    vec_QQ c2 = c;
    mat_ZZ num2 = num;

    for (int i = 0; i < c2.length(); ++i) {
	c2[i].canonicalize();
	if (c2[i].n != 0)
	    continue;

	if (i < c2.length()-1) {
	    num2[i] = num2[c2.length()-1];
	    c2[i] = c2[c2.length()-1];
	}
	num2.SetDims(num2.NumRows()-1, num2.NumCols());
	c2.SetLength(c2.length()-1);
	--i;
    }
    reduce(c2, num2, den_f);
}

void reducer::reduce(const vec_QQ& c, const mat_ZZ& num, const mat_ZZ& den_f)
{
    assert(c.length() == num.NumRows());
    unsigned len = den_f.NumRows();  // number of factors in den
    vec_QQ c2 = c;

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
    c2 *= sign;

    int only_param = 0;	    // k-r-s from text
    int no_param = 0;	    // r from text
    for (int k = 0; k < len; ++k) {
	if (den_p[k] == 0)
	    ++no_param;
	else if (den_s[k] == 0)
	    ++only_param;
    }
    if (no_param == 0) {
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
	for (int i = 0; i < num_s.length(); ++i) {
	    zz2value(num_s[i], tz);
	    n[i] = new dpoly(no_param, tz);
	    /* Search for other numerator (j) with same num_p.
	     * If found, replace a[j]/b[j] * n[j] and a[i]/b[i] * n[i]
	     * by 1/(b[j]*b[i]/g) * (a[j]*b[i]/g * n[j] + a[i]*b[j]/g * n[i])
	     * where g = gcd(b[i], b[j].
	     */
	    for (int j = 0; j < i; ++j) {
		if (num_p[i] != num_p[j])
		    continue;
		ZZ g = GCD(c2[i].d, c2[j].d);
		zz2value(c2[j].n * c2[i].d/g, tz);
		*n[j] *= tz;
		zz2value(c2[i].n * c2[j].d/g, tz);
		*n[i] *= tz;
		*n[j] += *n[i];
		c2[j].n = 1;
		c2[j].d *= c2[i].d/g;
		delete n[i];
		if (i < num_s.length()-1) {
		    num_s[i] = num_s[num_s.length()-1];
		    num_p[i] = num_p[num_s.length()-1];
		    c2[i] = c2[num_s.length()-1];
		}
		num_s.SetLength(num_s.length()-1);
		c2.SetLength(c2.length()-1);
		num_p.SetDims(num_p.NumRows()-1, num_p.NumCols());
		--i;
		break;
	    }
	}
	zz2value(den_s[k], tz);
	dpoly D(no_param, tz, 1);
	for ( ; ++k < len; )
	    if (den_p[k] == 0) {
		zz2value(den_s[k], tz);
		dpoly fact(no_param, tz, 1);
		D *= fact;
	    }

	if (no_param + only_param == len) {
	    vec_QQ q;
	    q.SetLength(num_s.length());
	    for (int i = 0; i < num_s.length(); ++i) {
		mpq_set_si(tcount, 0, 1);
		n[i]->div(D, tcount, 1);

		value2zz(mpq_numref(tcount), q[i].n);
		value2zz(mpq_denref(tcount), q[i].d);
		q[i] *= c2[i];
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
	    value_set_si(tz, 0);
	    dpoly one(no_param, tz);
	    dpoly_r *r = NULL;

	    for (k = 0; k < len; ++k) {
		if (den_s[k] == 0 || den_p[k] == 0)
		    continue;

		zz2value(den_s[k], tz);
		dpoly pd(no_param-1, tz, 1);

		int l;
		for (l = 0; l < k; ++l)
		    if (den_r[l] == den_r[k])
			break;

		if (!r)
		    r = new dpoly_r(one, pd, l, len);
		else {
		    dpoly_r *nr = new dpoly_r(r, pd, l, len);
		    delete r;
		    r = nr;
		}
	    }

	    vec_QQ factor;
	    factor.SetLength(c2.length());
	    int common = pden.NumRows();
	    dpoly_r *rc = r->div(D);
	    for (int i = 0; i < num_s.length(); ++i) {
		factor[i].d = c2[i].d;
		factor[i].d *= rc->denom;
	    }

	    dpoly_r_scanner scanner(n, num_s.length(), rc, len);
	    int rows;
	    while (scanner.next()) {
		int i;
		for (i = 0; i < num_s.length(); ++i)
		    if (scanner.coeff[i] != 0)
			break;
		if (i == num_s.length())
		    continue;
		rows = common;
		pden.SetDims(rows, pden.NumCols());
		for (int k = 0; k < rc->dim; ++k) {
		    int n = scanner.powers[k];
		    if (n == 0)
			continue;
		    pden.SetDims(rows+n, pden.NumCols());
		    for (int l = 0; l < n; ++l)
			pden[rows+l] = den_r[k];
		    rows += n;
		}
		/* The denominators in factor are kept constant
		 * over all iterations of the enclosing while loop.
		 * The rational numbers in factor may therefore not be
		 * canonicalized.  Some may even be zero.
		 */
		for (int i = 0; i < num_s.length(); ++i) {
		    factor[i].n = c2[i].n;
		    factor[i].n *= scanner.coeff[i];
		}
		reduce_canonical(factor, num_p, pden);
	    }

	    delete rc;
	    delete r;
	}
	for (int i = 0; i < num_s.length(); ++i)
	    delete n[i];
    }
}

void reducer::handle(const mat_ZZ& den, Value *V, const QQ& c,
		     unsigned long det, barvinok_options *options)
{
    vec_QQ vc;

    Matrix *points = Matrix_Alloc(det, dim);
    Matrix* Rays = zz2matrix(den);
    lattice_points_fixed(V, V, Rays, Rays, points, det);
    Matrix_Free(Rays);
    matrix2zz(points, vertex, points->NbRows, points->NbColumns);
    Matrix_Free(points);

    vc.SetLength(vertex.NumRows());
    for (int i = 0; i < vc.length(); ++i)
	vc[i] = c;

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

void icounter::base(const QQ& c, const vec_ZZ& num, const mat_ZZ& den_f)
{
    int r;
    unsigned len = den_f.NumRows();  // number of factors in den
    vec_ZZ den_s;
    den_s.SetLength(len);
    assert(num.length() == 1);
    ZZ num_s = num[0];
    for (r = 0; r < len; ++r)
	den_s[r] = den_f[r][0];
    int sign = (len % 2) ? -1 : 1;

    zz2value(num_s, tz);
    dpoly n(len, tz);
    zz2value(den_s[0], tz);
    dpoly D(len, tz, 1);
    for (int k = 1; k < len; ++k) {
	zz2value(den_s[k], tz);
	dpoly fact(len, tz, 1);
	D *= fact;
    }
    mpq_set_si(tcount, 0, 1);
    n.div(D, tcount, 1);
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
    assert(num.length() == 1);
    ZZ num_s = num[0];

    for (r = 0; r < len; ++r)
	den_s[r] = den_f[r][0];
    int sign = (len % 2) ? -1 : 1;

    zz2value(num_s, tz);
    dpoly n(len, tz);
    zz2value(den_s[0], tz);
    dpoly D(len, tz, 1);
    for (int k = 1; k < len; ++k) {
	zz2value(den_s[k], tz);
	dpoly fact(len, tz, 1);
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
