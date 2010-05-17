#include <barvinok/options.h>
#include "barvinok_union_options.h"

struct isl_arg union_options_arg[] = {
ISL_ARG_CHILD(struct union_options, barvinok, NULL, barvinok_options_arg, NULL)
ISL_ARG_BOOL(struct union_options, series, 's', "series", 0,
	"compute rational generating function")
ISL_ARG_END
};

ISL_ARG_DEF(union_options, struct union_options, union_options_arg)
