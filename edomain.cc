#include "fdstream.h"
#include "edomain.h"
#include "evalue_util.h"

using std::endl;

void EDomain::print(FILE *out, char **p)
{
    fdostream os(dup(fileno(out)));
    for (int i = 0; i < floors.size(); ++i) {
	os << "floor " << i << ": [";
	evalue_print(os, floors[i], p);
	os << "]" << endl;
    }
    Polyhedron_Print(out, P_VALUE_FMT, D);
}
