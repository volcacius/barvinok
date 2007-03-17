#ifndef VERIFY_H
#define VERIFY_H

#include <barvinok/polylib.h>
#include "argp.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct verify_options {
    int verify;
    int print_all;
    int continue_on_error;
    int m;
    int M;

    /* "generated" options */
    int st;
    char **params;
    struct barvinok_options *barvinok;
};

extern struct argp verify_argp;
void verify_options_set_range(struct verify_options *options, int dim);

#if defined(__cplusplus)
}
#endif

#endif
