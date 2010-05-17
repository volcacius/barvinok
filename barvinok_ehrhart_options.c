#include <barvinok/options.h>
#include "barvinok_ehrhart_options.h"

struct isl_arg ehrhart_options_arg[] = {
ISL_ARG_CHILD(struct ehrhart_options, convert, NULL,
	convert_options_arg, "output conversion")
ISL_ARG_CHILD(struct ehrhart_options, barvinok, NULL, barvinok_options_arg, NULL)
ISL_ARG_BOOL(struct ehrhart_options, series, 's', "series", 0,
	"compute rational generating function")
ISL_ARG_END
};

ISL_ARG_DEF(ehrhart_options, struct ehrhart_options, ehrhart_options_arg)
