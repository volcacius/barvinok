#include <vector>
#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/evalue.h>

struct EDomain {
    Polyhedron		*D;
    Vector		*sample;
    std::vector<evalue *>	floors;

    EDomain(Polyhedron *D) {
	this->D = Polyhedron_Copy(D);
	sample = NULL;
    }
    EDomain(Polyhedron *D, std::vector<evalue *>floors) {
	this->D = Polyhedron_Copy(D);
	add_floors(floors);
	sample = NULL;
    }
    EDomain(EDomain *ED) {
	this->D = Polyhedron_Copy(ED->D);
	add_floors(ED->floors);
	sample = NULL;
    }
    EDomain(Polyhedron *D, EDomain *ED, std::vector<evalue *>floors) {
	this->D = Polyhedron_Copy(D);
	add_floors(ED->floors);
	add_floors(floors);
	sample = NULL;
    }
    void add_floors(std::vector<evalue *>floors) {
	for (int i = 0; i < floors.size(); ++i) {
	    evalue *f = new evalue;
	    value_init(f->d);
	    evalue_copy(f, floors[i]);
	    this->floors.push_back(f);
	}
    }
    int find_floor(evalue *needle) {
	for (int i = 0; i < floors.size(); ++i)
	    if (eequal(needle, floors[i]))
		return i;
	return -1;
    }
    void print(FILE *out, char **p);
    ~EDomain() {
	for (int i = 0; i < floors.size(); ++i) {
	    free_evalue_refs(floors[i]);
	    delete floors[i];
	}
	Polyhedron_Free(D);
	if (sample)
	    Vector_Free(sample);
    }
};
