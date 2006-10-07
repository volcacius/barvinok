#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/options.h>

struct decomposer {
    void decompose(Polyhedron *C, barvinok_options *options);
    virtual void handle(Polyhedron *P, int sign) = 0;
};

struct polar_decomposer : public decomposer {
    void decompose(Polyhedron *C, barvinok_options *options);
    virtual void handle(Polyhedron *P, int sign);
    virtual void handle_polar(Polyhedron *P, int sign) = 0;
};

struct vertex_decomposer {
    Polyhedron *P;
    unsigned nbV;	// number of vertices
    Param_Vertices *V;	// current vertex
    int vert;		// current vertex index
    polar_decomposer *pd;

    vertex_decomposer(Polyhedron *P, unsigned nbV, polar_decomposer *pd) : 
			P(P), nbV(nbV), pd(pd) {}
    void decompose_at_vertex(Param_Vertices *V, int _i, barvinok_options *options);
};

#endif
