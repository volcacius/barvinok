#include <isl_arg.h>
#include "evalue_convert.h"
#include "verify.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct options {
	struct convert_options   *convert;
	struct verify_options    *verify;
	char *var_list;
};

ISL_ARG_DECL(options, struct options, options_arg)

#if defined(__cplusplus)
}
#endif
