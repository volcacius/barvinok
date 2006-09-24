#include <iostream>
#include <vector>
#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/evalue.h>
#include <barvinok/options.h>

struct EDomain_floor {
    int	     refcount;
    evalue  *e;
    Vector  *v;

    EDomain_floor(const evalue *f, int dim);
    ~EDomain_floor() {
	free_evalue_refs(e);
	delete e;
	Vector_Free(v);
    }
    EDomain_floor *ref() {
	++refcount;
	return this;
    }
    static void unref(EDomain_floor* floor) {
	if (!--floor->refcount)
	    delete floor;
    }
    void print(std::ostream& os, char **p) const;
};

struct EDomain {
    Polyhedron		*D;
    Vector		*sample;
    std::vector<EDomain_floor *>	floors;

    EDomain(Polyhedron *D) {
	this->D = Polyhedron_Copy(D);
	sample = NULL;
    }
    EDomain(Polyhedron *D, std::vector<EDomain_floor *>floors) {
	this->D = Polyhedron_Copy(D);
	add_floors(floors);
	sample = NULL;
    }
    EDomain(EDomain *ED) {
	this->D = Polyhedron_Copy(ED->D);
	add_floors(ED->floors);
	sample = NULL;
    }
    EDomain(Polyhedron *D, EDomain *ED, std::vector<EDomain_floor *>floors) {
	this->D = Polyhedron_Copy(D);
	add_floors(ED->floors);
	add_floors(floors);
	sample = NULL;
    }
    void add_floors(std::vector<EDomain_floor *>floors) {
	for (int i = 0; i < floors.size(); ++i)
	    this->floors.push_back(floors[i]->ref());
    }
    int find_floor(evalue *needle) const {
	for (int i = 0; i < floors.size(); ++i)
	    if (eequal(needle, floors[i]->e))
		return i;
	return -1;
    }
    void print(FILE *out, char **p);
    void print_constraints(std::ostream& os, char **p,
			   barvinok_options *options) const;
    ~EDomain() {
	for (int i = 0; i < floors.size(); ++i)
	    EDomain_floor::unref(floors[i]);
	Polyhedron_Free(D);
	if (sample)
	    Vector_Free(sample);
    }

    unsigned dimension() const {
	return D->Dimension - floors.size();
    }

    Matrix *add_ge_constraint(evalue *constraint,
				 std::vector<EDomain_floor *>& new_floors) const;
};

int evalue2constraint(EDomain *D, evalue *E, Value *cons, int len);
