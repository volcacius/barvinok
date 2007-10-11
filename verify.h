#ifndef VERIFY_H
#define VERIFY_H

#include <barvinok/polylib.h>
#include "argp.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct verify_options {
    int verify;
    int exact;
    int print_all;
    int continue_on_error;
    int m;
    int M;

    /* "generated" options */
    int st;
    const char *const *params;
    struct barvinok_options *barvinok;
};

extern struct argp verify_argp;
void verify_options_set_range(struct verify_options *options, int dim);

Polyhedron *check_poly_context_scan(Polyhedron *P, Polyhedron **C,
				    unsigned nparam,
				    const struct verify_options *options);
void check_poly_init(Polyhedron *C, struct verify_options *options);
void check_poly_print(int ok, int nparam, Value *z,
		      Value want_n, Value want_d,
		      Value got_n, Value got_d,
		      const char *op, const char *check,
		      const char *long_op,
		      const struct verify_options *options);

struct check_poly_data;
typedef int (*check_poly_fun)(const struct check_poly_data *data,
			      int nparam, Value *z,
		 	      const struct verify_options *options);
struct check_poly_data {
    Value	    *z;
    check_poly_fun   check;
};

int check_poly(Polyhedron *CS, const struct check_poly_data *data,
	       int nparam, int pos, Value *z,
	       const struct verify_options *options);

#if defined(__cplusplus)
}
#endif

#endif
