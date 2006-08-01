#include "genfun_constructor.h"

gf_base *gf_base::create(Polyhedron *context, unsigned dim, unsigned nparam)
{
    gf_base *red;
#ifdef USE_INCREMENTAL_BF
    red = new partial_bfcounter(context, dim, nparam);
#elif defined USE_INCREMENTAL_DF
    red = new partial_ireducer(context, dim, nparam);
#else
    red = new partial_reducer(context, dim, nparam);
#endif
    return red;
}

void partial_ireducer::base(QQ& c, const vec_ZZ& num, const mat_ZZ& den_f)
{
    gf->add(c, num, den_f);
}

void partial_reducer::split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
			    mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r)
{
    unsigned len = den_f.NumRows();  // number of factors in den
    unsigned nvar = tmp.length();

    den_s.SetLength(len);
    den_r.SetDims(len, lower);

    for (int r = 0; r < len; ++r) {
	for (int k = 0; k < nvar; ++k)
	    tmp[k] = den_f[r][k];
	den_s[r] = tmp * lambda;

	for (int k = nvar; k < dim; ++k)
	    den_r[r][k-nvar] = den_f[r][k];
    }

    for (int k = 0; k < nvar; ++k)
	tmp[k] = num[k];
    num_s = tmp *lambda;
    num_p.SetLength(lower);
    for (int k = nvar ; k < dim; ++k)
	num_p[k-nvar] = num[k];
}

void partial_reducer::base(QQ& c, const vec_ZZ& num, const mat_ZZ& den_f)
{
    gf->add(c, num, den_f);
}

void partial_bfcounter::base(mat_ZZ& factors, bfc_vec& v)
{
    mat_ZZ den;
    unsigned nf = factors.NumRows();

    for (int i = 0; i < v.size(); ++i) {
	bfc_term* bfct = static_cast<bfc_term *>(v[i]);
	den.SetDims(0, lower);
	int total_power = 0;
	int p = 0;
	for (int j = 0; j < nf; ++j) {
	    total_power += v[i]->powers[j];
	    den.SetDims(total_power, lower);
	    for (int k = 0; k < v[i]->powers[j]; ++k)
		den[p++] = factors[j];
	}
	for (int j = 0; j < v[i]->terms.NumRows(); ++j)
	    gf->add(bfct->c[j], v[i]->terms[j], den);
	delete v[i];
    }
}
