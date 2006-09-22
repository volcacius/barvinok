#include <vector>
#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/evalue.h>

struct EDomain_floor {
    int	     refcount;
    evalue  *e;

    EDomain_floor(const evalue *f) {
	e = new evalue;
	value_init(e->d);
	evalue_copy(e, f);
	refcount = 1;
    }
    ~EDomain_floor() {
	free_evalue_refs(e);
	delete e;
    }
    EDomain_floor *ref() {
	++refcount;
	return this;
    }
    static void unref(EDomain_floor* floor) {
	if (!--floor->refcount)
	    delete floor;
    }
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
    int find_floor(evalue *needle) {
	for (int i = 0; i < floors.size(); ++i)
	    if (eequal(needle, floors[i]->e))
		return i;
	return -1;
    }
    void print(FILE *out, char **p);
    ~EDomain() {
	for (int i = 0; i < floors.size(); ++i)
	    EDomain_floor::unref(floors[i]);
	Polyhedron_Free(D);
	if (sample)
	    Vector_Free(sample);
    }
};
