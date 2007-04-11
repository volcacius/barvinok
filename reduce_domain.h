#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_options;

Vector *inner_point(Polyhedron *P);
int is_internal(Vector *point, Value *constraint);
Polyhedron *reduce_domain(Polyhedron *D, Matrix *CT, Polyhedron *CEq,
			  Polyhedron **fVD, int nd,
			  struct barvinok_options *options);

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define FORALL_REDUCED_DOMAIN(PP,CT,CEq,nd,options,i,D,rVD)		    \
	do {								    \
	    Param_Domain *D;						    \
	    Polyhedron *rVD;						    \
	    int i;							    \
	    int _frd_i;							    \
	    Polyhedron **_frd_fVD;					    \
	    Param_Domain *_frd_next;					    \
	    if (nd < 0)							    \
		for (nd = 0, D = PP->D; D; ++nd, D = D->next);		    \
	    _frd_fVD = ALLOCN(Polyhedron *, nd);			    \
	    for (i = 0, D = PP->D; D; D = _frd_next) {			    \
		rVD = reduce_domain(D->Domain, CT, CEq, _frd_fVD,	    \
				    i, options);			    \
		_frd_next = D->next;					    \
		if (!rVD)						    \
		    continue;						    \
		{
#define FORALL_REDUCED_DOMAIN_RESET(nd)					    \
	do {								    \
	    int _frd_i;							    \
	    for (_frd_i = 0; _frd_i < (nd); ++_frd_i)			    \
		Domain_Free(_frd_fVD[_frd_i]);				    \
	    free(_frd_fVD);						    \
	} while(0);
#define END_FORALL_REDUCED_DOMAIN					    \
		}							    \
		++i;							    \
	    }								    \
	    FORALL_REDUCED_DOMAIN_RESET(i)				    \
	    nd = i;							    \
	} while(0);


#if defined(__cplusplus)
}
#endif
