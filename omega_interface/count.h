#include <barvinok/barvinok.h>
#include <omega.h>

#define COUNT_RELATION_BARVINOK	0
#define COUNT_RELATION_PARKER	1

evalue *count_relation(Relation& r, int method);
evalue *rank_relation(Relation& r);
evalue *count_lexsmaller(Relation& r, Relation& domain);

evalue *barvinok_enumerate_parker(Polyhedron *D,
					unsigned nvar, unsigned nparam,
					unsigned MaxRays);
