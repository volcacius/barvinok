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

struct piecewise_lst {
    const GiNaC::exvector vars;
    std::vector<guarded_lst> list;
    /*  0: just collect terms
     *  1: remove obviously smaller terms (maximize)
     * -1: remove obviously bigger terms (minimize)
     */
    int sign;

    piecewise_lst(const GiNaC::exvector& vars, int sign = 0) :
			vars(vars), sign(sign) {}
    void add_guarded_lst(Polyhedron *D, GiNaC::lst coeffs);
    piecewise_lst& combine(const piecewise_lst& other);
    void maximize();
    void minimize();
    void simplify_domains(Polyhedron *ctx, unsigned MaxRays);
    GiNaC::numeric evaluate(const GiNaC::exvector& values);
    void add(const GiNaC::ex& poly);

    ~piecewise_lst() {
	free_list_domains();
    }
private:
    /* We don't allow making copies (yet), because we would have
     * to make copies of all polyhedra.
     */
    piecewise_lst(const piecewise_lst&) {}
    void free_list_domains() {
	for (int i = 0; i < list.size(); ++i)
	    Domain_Free(list[i].first);
    }
};

std::ostream & operator<< (std::ostream & os, const piecewise_lst & pl);

}

#endif
