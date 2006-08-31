#ifndef BERNSTEIN_PIECWISE_LST_H
#define BERNSTEIN_PIECWISE_LST_H

#include <iostream>
#include <vector>
#include <gmp.h>
#include <ginac/ginac.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

namespace bernstein {

typedef std::pair< Polyhedron *, GiNaC::lst > guarded_lst;

struct piecewise_lst_s {
    const GiNaC::exvector vars;
    std::vector<guarded_lst> list;

    piecewise_lst_s::piecewise_lst_s(const GiNaC::exvector& vars) : vars(vars) {}
    piecewise_lst_s& combine(const piecewise_lst_s& other);
    void maximize();
    void simplify_domains(Polyhedron *ctx);
    GiNaC::numeric evaluate(const GiNaC::exvector& values);
    void add(const GiNaC::ex& poly);

    ~piecewise_lst_s() {
	free_list_domains();
    }
private:
    /* We don't allow making copies (yet), because we would have
     * to make copies of all polyhedra.
     */
    piecewise_lst_s::piecewise_lst_s(const piecewise_lst_s&) {}
    void free_list_domains() {
	for (int i = 0; i < list.size(); ++i)
	    Domain_Free(list[i].first);
    }
};

std::ostream & operator<< (std::ostream & os, const piecewise_lst_s & pl);

// typedef GiNaC::structure<piecewise_lst_s> piecewise_lst;

}

#endif
