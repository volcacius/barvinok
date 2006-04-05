#include <bernstein/bernstein.h>
#include <bernstein/maximize.h>
#include <bernstein/piecewise_lst.h>

#define MAXRAYS    POL_NO_DUAL

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

std::ostream & operator<< (std::ostream & os, const piecewise_lst_s & pl)
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
					const std::vector<guarded_lst>& other)
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
	    lst list = one[i].second;
	    lst::const_iterator k;
	    for (k = other[j].second.begin(); k != other[j].second.end(); ++k)
		list.append(*k);
	    //std::copy(other[j].second.begin(), other[j].second.end(), 
		      //std::back_insert_iterator<lst>(list));
	    newlist.push_back(guarded_lst(d, list.sort().unique()));
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

piecewise_lst_s& piecewise_lst_s::combine(const piecewise_lst_s& other)
{
    assert(vars == other.vars);
    list = bernstein::combine(list, other.list);
    return *this;
}


void piecewise_lst_s::maximize()
{
    exvector params;
    for (int i = 0; i < list.size(); ++i) {
	list[i].second = bernstein::maximize(list[i].first, list[i].second, vars);
    }
}

void piecewise_lst_s::simplify_domains(Polyhedron *ctx)
{
    for (int i = 0; i < list.size(); ++i) {
	Polyhedron *D = list[i].first;
	list[i].first = DomainSimplify(D, ctx, MAXRAYS);
	Domain_Free(D);
    }
}

}
