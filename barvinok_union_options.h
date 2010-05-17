#include <isl_arg.h>
#include <barvinok/options.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct union_options {
	int series;
	struct barvinok_options *barvinok;
};

ISL_ARG_DECL(union_options, struct union_options, union_options_arg)

#if defined(__cplusplus)
}
#endif
