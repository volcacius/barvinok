#include <isl_arg.h>
#include "evalue_convert.h"
#include "verify.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct options {
	struct verify_options    *verify;
	long split;
	int lower;
	long iterate;
};

ISL_ARG_DECL(options, struct options, options_arg)

#if defined(__cplusplus)
}
#endif
