#include <assert.h>
#include <cln/cln.h>

#include <bernstein/bernstein.h>
#include <bernstein/maximize.h>
#include <bernstein/piecewise_lst.h>

using std::cerr;
using std::endl;

using namespace GiNaC;

namespace bernstein {

static void print_max(std::ostream & os, lst l)
{
    for (int i = 1; i < l.nops(); ++i)
	os << "max(";
    if (l.nops() > 0)
	l.op(0).print(print_csrc(os), 0);
    for (int i = 1; i < l.nops(); ++i) {
	os << ",";
	l.op(i).print(print_csrc(os), 0);
	os << ")";
    }
}

static void printpoly(std::ostream& o, Polyhedron *D, const exvector& p)
{
    for (int i = 0; i < D->NbConstraints; ++i) {
	int first = 1;
	if (i)
	    o << " && ";
	for (int j = 0; j < D->Dimension; ++j) {
	    if (value_zero_p(D->Constraint[i][1+j]))
		continue;
	    if (!first)
		o << " + ";
	    o << VALUE_TO_INT(D->Constraint[i][1+j]) << "*" << p[j];
	    first = 0;
	}
	if (value_notzero_p(D->Constraint[i][1+D->Dimension])) {
	    if (!first)
		o << " + ";
	    o << VALUE_TO_INT(D->Constraint[i][1+D->Dimension]);
	}
	if (value_zero_p(D->Constraint[i][0]))
	    o << " == 0";
	else
	    o << " >= 0";
    }
}

static void printdomain(std::ostream& o, Polyhedron *D, const exvector& p)
{
    if (!D->next)
	printpoly(o, D, p);
    else for (Polyhedron *P = D; P; P = P->next) {
	o << "(";
	printpoly(o, P, p);
	o << ")";
	if (P->next)
	    o << " || ";
    }
}

std::ostream & operator<< (std::ostream & os, const piecewise_lst & pl)
{
    if (pl.list.size() == 1 && universeQ(pl.list[0].first))
	print_max(os, pl.list[0].second);
    else {
	for (int i = 0; i < pl.list.size(); ++i) {
	    os << "(";
	    printdomain(os, pl.list[i].first, pl.vars);
	    os << ") ? (";
	    print_max(os, pl.list[i].second);
	    os << ") : ";
	}
	os << "0";
    }
    return os;
}

static std::vector<guarded_lst> combine(const std::vector<guarded_lst>& one,
					const std::vector<guarded_lst>& other,
					const GiNaC::exvector& vars, int sign)
{
    Polyhedron *fd;
    std::vector<guarded_lst> newlist;
    for (int j = 0; j < other.size(); ++j) {
	assert(one.size() >= 1);
	fd = DomainDifference(other[j].first, one[0].first, 0);
	if (!emptyQ(fd))
	    for (int i = 1; i < one.size(); ++i) {
		Polyhedron *t = fd;
		fd = DomainDifference(fd, one[i].first, 0);
		Domain_Free(t);
		if (emptyQ(fd))
		    break;
	    }
	if (emptyQ(fd)) {
	    Domain_Free(fd);
	    continue;
	}
	newlist.push_back(guarded_lst(fd, other[j].second));
    }
    for (int i = 0; i < one.size(); ++i) {
	fd = one[i].first;
	for (int j = 0; j < other.size(); ++j) {
	    Polyhedron *t, *d;
	    d = DomainIntersection(other[j].first, one[i].first, 0);
	    if (emptyQ(d)) {
		Domain_Free(d);
		continue;
	    }
	    t = fd;
	    fd = DomainDifference(fd, other[j].first, 0);
	    if (t != one[i].first)
		Domain_Free(t);
	    lst list = remove_redundants(d, one[i].second, other[j].second,
					 vars, sign);
	    newlist.push_back(guarded_lst(d, list));
	}
	if (!emptyQ(fd))
	    newlist.push_back(guarded_lst(fd, one[i].second));
	else
	    Domain_Free(fd);
	if (fd != one[i].first)
	    Domain_Free(one[i].first);
    }
    return newlist;
}
 
using std::ostream;

ostream & operator<< (ostream & os, const exvector & v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
	if (i)
	    os << ", ";
	os << v[i];
    }
    os << "]";
    return os;
}

void piecewise_lst::add_guarded_lst(Polyhedron *D, GiNaC::lst coeffs)
{
    coeffs = remove_redundants(D, coeffs, coeffs, vars, sign);
    assert(coeffs.nops() > 0);
    list.push_back(guarded_lst(D, coeffs));
}

piecewise_lst& piecewise_lst::combine(const piecewise_lst& other)
{
    assert(vars == other.vars);
    assert(sign == other.sign);
    list = bernstein::combine(list, other.list, vars, sign);
    return *this;
}


void piecewise_lst::maximize()
{
    exvector params;
    for (int i = 0; i < list.size(); ++i) {
	list[i].second = bernstein::maximize(list[i].first, list[i].second, vars);
    }
}

void piecewise_lst::minimize()
{
    exvector params;
    for (int i = 0; i < list.size(); ++i) {
	list[i].second = bernstein::minimize(list[i].first, list[i].second, vars);
    }
}

void piecewise_lst::simplify_domains(Polyhedron *ctx, unsigned MaxRays)
{
    for (int i = 0; i < list.size(); ++i) {
	Polyhedron *D = list[i].first;
	list[i].first = DomainSimplify(D, ctx, MaxRays);
	Domain_Free(D);
    }
}

numeric piecewise_lst::evaluate(const exvector& values)
{
    Value *v = new Value[values.size()];
    numeric result;

    for (int i = 0; i < values.size(); ++i) {
	value_init(v[i]);
	assert(is_a<numeric>(values[i]));
	numeric2value(ex_to<numeric>(values[i]), v[i]);
    }
    result = evaluate(values, values.size(), v);
    for (int i = 0; i < values.size(); ++i)
	value_clear(v[i]);
    delete [] v;
    return result;
}

void piecewise_lst::evaluate(int n, Value *v, Value *num, Value *den)
{
    exvector values(n);
    numeric result;
    for (int i = 0; i < values.size(); ++i)
	values[i] = value2numeric(v[i]);
    result = evaluate(values, values.size(), v);
    numeric2value(result.numer(), *num);
    numeric2value(result.denom(), *den);
}

numeric piecewise_lst::evaluate(const exvector& values, int n, Value *v)
{
    numeric result = 0;

    for (int i = 0; i < list.size(); ++i) {
	if (!in_domain(list[i].first, v))
	    continue;
	exmap m;
	for (int j = 0; j < n; ++j)
	    m[vars[j]] = values[j];
	ex ex_val = list[i].second.subs(m);
	assert(is_a<lst>(ex_val));
	lst val = ex_to<lst>(ex_val);;
	ex max = val.op(0);
	for (int j = 1; j < val.nops(); ++j)
	    if (val.op(j) > max)
		max = val.op(j);
	assert(is_a<numeric>(max));
	result = ex_to<numeric>(max);
	break;
    }
    return result;
}

void piecewise_lst::add(const GiNaC::ex& poly)
{
    for (int i = 0; i < list.size(); ++i) {
	lst::const_iterator j;
	lst new_coeff;
	for (j = list[i].second.begin(); j != list[i].second.end(); ++j)
	    new_coeff.append(*j + poly);
	list[i].second = new_coeff;
    }
}

int piecewise_lst::is_equal(const piecewise_lst& other) const
{
    if (list.size() != other.list.size())
	return 0;

    for (int i = 0; i < list.size(); ++i) {
	int j;
	lst l1, l2;
	for (j = 0; j < other.list.size(); ++j) {
	    if (!PolyhedronIncludes(list[i].first, other.list[j].first))
		continue;
	    if (!PolyhedronIncludes(other.list[j].first, list[i].first))
		continue;
	    break;
	}
	if (j >= other.list.size())
	    return 0;
	l1 = list[i].second;
	l2 = other.list[j].second;
	l1.sort();
	l2.sort();
	if (!l1.is_equal(l2))
	    return 0;
    }

    return 1;
}

}
