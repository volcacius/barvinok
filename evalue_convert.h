#ifndef EVALUE_CONVERT
#define EVALUE_CONVERT

#include <barvinok/evalue.h>
#include "argp.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct convert_options {
    int range;
    int convert;
    int combine;
    int floor;
    int list;
    int latex;
};

int evalue_convert(evalue *EP, struct convert_options *options,
		   int verbose, unsigned nparam, char **params);

extern struct argp convert_argp;

#if defined(__cplusplus)
}
#endif

#endif
