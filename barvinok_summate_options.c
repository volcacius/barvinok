#include "barvinok_summate_options.h"

struct isl_arg options_arg[] = {
ISL_ARG_CHILD(struct options, verify, NULL,
	verify_options_arg, "verification")
ISL_ARG_CHILD(struct options, convert, NULL,
	convert_options_arg, "output conversion")
ISL_ARG_STR(struct options, var_list, 0, "variables", "list", NULL,
	"comma separated list of variables over which to sum")
ISL_ARG_END
};

ISL_ARG_DEF(options, struct options, options_arg)
