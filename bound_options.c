#include "bound_options.h"

struct isl_arg options_arg[] = {
ISL_ARG_CHILD(struct options, verify, NULL,
	verify_options_arg, "verification")
ISL_ARG_CHILD(struct options, convert, NULL,
	convert_options_arg, "input conversion")
ISL_ARG_STR(struct options, var_list, 0, "variables", "list", NULL,
	"comma separated list of variables over which to sum")
ISL_ARG_LONG(struct options, split, 0, "split", 0, NULL)
ISL_ARG_OPT_LONG(struct options, iterate, 0, "iterate", 0, -1,
	"exact result by iterating over domain (of specified maximal size)")
ISL_ARG_BOOL(struct options, lower, 0, "lower", 0,
	"compute lower bound instead of upper bound")
ISL_ARG_END
};

ISL_ARG_DEF(options, struct options, options_arg)
