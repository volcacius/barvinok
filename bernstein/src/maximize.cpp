#include <bernstein/bernstein.h>
#include <bernstein/maximize.h>

using namespace GiNaC;

namespace bernstein {

static void domainVerticesAndRays(Polyhedron *P, GiNaC::matrix& VM, GiNaC::matrix& RM)
{
    POL_ENSURE_VERTICES(P);
    int nv = 0, nr = 0;
    for (int i = 0; i < P->NbRays; ++i) {
	if (value_zero_p(P->Ray[i][0]))
	    nr += 2;
	else if (value_zero_p(P->Ray[i][1+P->Dimension]))
	    ++nr;
	else
	    ++nv;
    }
    VM = GiNaC::matrix(nv, P->Dimension);
    RM = GiNaC::matrix(nr, P->Dimension);
    nv = nr = 0;
    for (int i = 0; i < P->NbRays; ++i) {
	if (value_zero_p(P->Ray[i][1+P->Dimension])) {
	    for (int j = 0; j < P->Dimension; ++j)
		RM(nr,j) = value2numeric(P->Ray[i][1+j]);
	    ++nr;
	    if (value_notzero_p(P->Ray[i][0]))
		continue;
	    for (int j = 0; j < P->Dimension; ++j)
		RM(nr,j) = -RM(nr-1,j);
	    ++nr;
	    continue;
	}
	for (int j = 0; j < P->Dimension; ++j)
	    VM(nv,j) = value2numeric(P->Ray[i][1+j]) / 
					value2numeric(P->Ray[i][1+P->Dimension]);
	++nv;
    }
}

static void find_lst_minmax(const lst& l, ex& m, ex& M)
{
    lst::const_iterator i = l.begin();
    m = M = *i;
    for (++i; i != l.end(); ++i) {
	if (*i > M)
	    M = *i;
	else if (*i < m)
	    m = *i;
    }
}

static bool is_linear_constraint(ex constraint, const exvector& vars, matrix& linear)
{
    linear = matrix(vars.size()+1, 1);
    ex cst = constraint;
    for (int i = 0; i < vars.size(); ++i) {
	if (cst.degree(vars[i]) > 1)
	    return false;
	ex c = cst.coeff(vars[i], 1);
	if (!is_a<numeric>(c))
	    return false;
	linear(i,0) = c;
	cst = cst.coeff(vars[i], 0);
    }
    if (!is_a<numeric>(cst))
	return false;
    linear(vars.size(),0) = cst;
    return true;
}

static void find_linear_minmax(Polyhedron *D, matrix constraint, ex& m, ex& M)
{
    // we only care about the sign anyway
    M = -1; m = 1;
    for (Polyhedron *P = D; P; P = P->next) {
	POL_ENSURE_VERTICES(P);
	for (int i = 0; i < P->NbRays; ++i) {
	    ex sum = 0;
	    assert(value_one_p(P->Ray[i][0]));
	    for (int j = 0; j < P->Dimension+1; ++j)
		sum += value2numeric(P->Ray[i][1+j]) * constraint(j,0);
	    if (sum > M)
		M = sum;
	    if (sum < m)
		m = sum;
	}
    }
}

enum poly_sign {
    zero,
    positive,
    negative,
    unknown,
};

static poly_sign polynomial_sign(ex poly, Polyhedron *D, const GiNaC::matrix& VM,
				 const GiNaC::matrix& RM, const exvector& vars);

/* bernstein requires a bounded domain, so if the domain is unbounded,
 * we need a different trick.
 * If a parameter N has a lowerbound l, then we write the polynomial as
 * q(N,M) * (N-l) + r(M)
 * if the signs of q and r can be determined and if they are equal,
 * then the whole polynomial has the same sign.
 */
static poly_sign polynomial_sign_with_rays(ex poly, Polyhedron *D, 
				 const GiNaC::matrix& VM,
				 const GiNaC::matrix& RM, const exvector& vars)
{
    int i;
    for (i = 0; i < vars.size(); ++i)
	if (poly.degree(vars[i]) > 0)
	    break;
    assert(i < vars.size());
    ex min = VM(0,i);
    for (int j = 1; j < VM.rows(); ++j)
	if (VM(j,i) < min)
	    min = VM(j,i);
    if (min < 0)
	return unknown;
    for (int j = 0; j < RM.rows(); ++j)
	if (RM(j,i) < 0)
	    return unknown;
    int d = poly.degree(vars[i]);
    ex cum = poly.coeff(vars[i],d);
    ex q = cum;

    for (--d; d >= 1; --d) {
	cum *= min;
	cum += poly.coeff(vars[i],d);
	q *= vars[i];
	q += cum;
    }
    ex r = cum * min + poly.coeff(vars[i], 0);
    poly_sign s1 = polynomial_sign(q.expand(), D, VM, RM, vars);
    if (s1 == unknown)
	return unknown;
    poly_sign s2 = polynomial_sign(r, D, VM, RM, vars);
    if (s2 == unknown)
	return unknown;
    if (s1 == zero)
	return s2;
    if (s2 == zero)
	return s1;
    if (s1 != s2)
	return unknown;
    return s1;
}

static ex substitute_equalities(ex poly, Polyhedron *D, const exvector& vars)
{
    if (D->NbEq == 0)
	return poly;
    lst replacements;
    for (int i = 0; i < D->NbEq; ++i) {
	int pos = First_Non_Zero(D->Constraint[i]+1, D->Dimension);
	assert(pos != -1);
	ex den = -value2numeric(D->Constraint[i][1+pos]);
	ex s = value2numeric(D->Constraint[i][1+D->Dimension]) / den;
	for (int j = pos+1; j < D->Dimension; ++j) {
	    if (value_zero_p(D->Constraint[i][1+j]))
		continue;
	    s += value2numeric(D->Constraint[i][1+j]) * vars[j] / den;
	}
	replacements.append(vars[pos] == s);
    }
    poly = poly.subs(replacements).expand();
    return poly;
}

static poly_sign polynomial_sign(ex poly, Polyhedron *D, const GiNaC::matrix& VM,
				 const GiNaC::matrix& RM, const exvector& vars)
{
    exvector params;
    ex minc, maxc;
    matrix linear;
    poly = substitute_equalities(poly, D, vars);
    if (is_a<numeric>(poly))
	minc = maxc = poly;
    else if (is_linear_constraint(poly, vars, linear))
	find_linear_minmax(D, linear, minc, maxc);
    else if (RM.rows() == 0) {
	lst coeffs = bernsteinExpansion(VM, poly, vars, params);
	find_lst_minmax(coeffs, minc, maxc);
    } else
	return polynomial_sign_with_rays(poly, D, VM, RM, vars);
    if (maxc <= 0 && minc >= 0)
	return zero;
    if (maxc <= 0)
	return negative;
    if (minc >= 0)
	return positive;
    return unknown;
}

GiNaC::lst maximize(Polyhedron *domain, GiNaC::lst coeffs,
		    const GiNaC::exvector& vars)
{
    GiNaC::matrix VM, RM;
    domainVerticesAndRays(domain, VM, RM);
    lst::const_iterator j, k;
    lst newlist;
    lst removed;
    for (j = coeffs.begin(); j != coeffs.end(); ++j) {
	if (find(removed.begin(), removed.end(), *j) != removed.end())
	    continue;
	bool needed = true;
	k = j; ++k;
	for (; k != coeffs.end(); ++k) {
	    if (find(removed.begin(), removed.end(), *k) != removed.end())
		continue;
	    ex diff = *j - *k;
	    poly_sign s = polynomial_sign(diff, domain, VM, RM, vars);
	    if (s == zero || s == negative)
		needed = false;
	    else if (s == positive)
		removed.append(*k);
	    if (!needed)
		break;
	}
	if (needed)
	    newlist.append(*j);
    }
    return newlist;
}

GiNaC::lst minimize(Polyhedron *domain, GiNaC::lst coeffs,
		    const GiNaC::exvector& vars)
{
    lst negcoeffs, negresult, result;
    lst::const_iterator j;
    for (j = coeffs.begin(); j != coeffs.end(); ++j)
	negcoeffs.append(-*j);
    negresult = maximize(domain, negcoeffs, vars);
    for (j = negresult.begin(); j != negresult.end(); ++j)
	result.append(-*j);
    return result;
}

}
