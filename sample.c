#include <assert.h>
#include <isl_set_polylib.h>
#include <barvinok/sample.h>

Vector *Polyhedron_Sample(Polyhedron *P, struct barvinok_options *options)
{
	int i;
	isl_ctx *ctx = isl_ctx_alloc();
	isl_space *dim;
	int nvar = P->Dimension;
	isl_basic_set *bset;
	isl_point *pnt;
	Vector *sample = NULL;

	dim = isl_space_set_alloc(ctx, 0, nvar);
	bset = isl_basic_set_new_from_polylib(P, dim);
	pnt = isl_basic_set_sample_point(bset);

	if (!isl_point_is_void(pnt)) {
		isl_int v;

		isl_int_init(v);
		sample = Vector_Alloc(1 + nvar);
		assert(sample);
		for (i = 0; i < nvar; ++i) {
			isl_point_get_coordinate(pnt, isl_dim_set, i, &v);
			isl_int_get_gmp(v, sample->p[i]);
		}
		value_set_si(sample->p[nvar], 1);
		isl_int_clear(v);
	}

	isl_point_free(pnt);

	isl_ctx_free(ctx);

	return sample;
}
