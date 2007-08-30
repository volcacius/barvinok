#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_options;

Polyhedron *true_context(Polyhedron *P, Polyhedron *C, unsigned MaxRays);
Vector *inner_point(Polyhedron *P);
int is_internal(Vector *point, Value *constraint);

Polyhedron *reduce_domain(Polyhedron *D, int nd,
			  Vector *inner, struct barvinok_options *options);

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define FORALL_REDUCED_DOMAIN(PP,C,nd,options,i,PD,rVD)			    \
	do {								    \
	    Param_Domain *PD;						    \
	    Polyhedron *rVD;						    \
	    int i;							    \
	    Param_Domain *_frd_next;					    \
	    Vector *_frd_inner = inner_point(C);			    \
	    if (nd < 0)							    \
		for (nd = 0, PD = PP->D; PD; ++nd, PD = PD->next);	    \
	    for (i = 0, PD = PP->D; PD; PD = _frd_next) {		    \
		rVD = reduce_domain(PD->Domain, nd, _frd_inner, options);   \
		_frd_next = PD->next;					    \
		if (!rVD)						    \
		    continue;						    \
		{
#define FORALL_REDUCED_DOMAIN_RESET					    \
	do {								    \
	    Vector_Free(_frd_inner);					    \
	} while(0);
#define END_FORALL_REDUCED_DOMAIN					    \
		}							    \
		++i;							    \
	    }								    \
	    FORALL_REDUCED_DOMAIN_RESET				    	    \
	    nd = i;							    \
	} while(0);


#if defined(__cplusplus)
}
#endif
