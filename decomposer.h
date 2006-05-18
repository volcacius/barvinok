#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

struct decomposer {
    void decompose(Polyhedron *C);
    virtual void handle(Polyhedron *P, int sign) = 0;
};

struct polar_decomposer : public decomposer {
    void decompose(Polyhedron *C, unsigned MaxRays);
    virtual void handle(Polyhedron *P, int sign);
    virtual void handle_polar(Polyhedron *P, int sign) = 0;
};

#endif
