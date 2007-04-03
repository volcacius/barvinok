#include <barvinok/polylib.h>

struct barvinok_options;

#if defined(__cplusplus)
extern "C" {
#endif

struct scale_data {
    Value   	det;
    int	    	save_approximation;
};

Polyhedron *Polyhedron_Flate(Polyhedron *P, unsigned nparam, int inflate,
			     unsigned MaxRays);

Polyhedron *scale_init(Polyhedron *P, Polyhedron *C, struct scale_data *scaling,
		       struct barvinok_options *options);
Polyhedron *scale(Param_Polyhedron *PP, Polyhedron *P,
		  struct scale_data *scaling, int free_P,
		  struct barvinok_options *options);
void scale_finish(evalue *e, struct scale_data *scaling,
		  struct barvinok_options *options);

#if defined(__cplusplus)
}
#endif
