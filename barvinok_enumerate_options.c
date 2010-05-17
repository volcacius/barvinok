#include "barvinok_enumerate_options.h"

struct isl_arg enumerate_options_arg[] = {
ISL_ARG_CHILD(struct enumerate_options, verify, NULL,
	verify_options_arg, "verification")
ISL_ARG_CHILD(struct enumerate_options, convert, NULL,
	convert_options_arg, "output conversion")
ISL_ARG_BOOL(struct enumerate_options, size, 'S', "size", 0, NULL)
ISL_ARG_BOOL(struct enumerate_options, series, 's', "series", 0,
	"compute rational generating function")
ISL_ARG_BOOL(struct enumerate_options, function, 'e', "explicit", 0,
	"convert rgf to psp")
ISL_ARG_END
};

ISL_ARG_DEF(enumerate_options, struct enumerate_options, enumerate_options_arg)
