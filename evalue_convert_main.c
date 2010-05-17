#include <barvinok/options.h>
#include "evalue_convert.h"
#include "evalue_read.h"

struct options {
	struct convert_options    *convert;
	struct barvinok_options   *barvinok;
};

struct isl_arg options_arg[] = {
ISL_ARG_CHILD(struct options, convert, NULL, convert_options_arg, NULL)
ISL_ARG_CHILD(struct options, barvinok, NULL, barvinok_options_arg, NULL)
ISL_ARG_END
};

ISL_ARG_DEF(options, struct options, options_arg)

int main(int argc, char **argv)
{
    evalue *EP;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct options *options = options_new_with_defaults();
    int printed;

    argc = options_parse(options, argc, argv, ISL_ARG_ALL);

    EP = evalue_read_from_file(stdin, NULL, &all_vars,
			       &nvar, &nparam, options->barvinok->MaxRays);
    assert(EP);

    printed =
	evalue_convert(EP, options->convert, options->barvinok->verbose,
			nvar+nparam, all_vars);

    if (!printed)
	print_evalue(stdout, EP, all_vars);

    evalue_free(EP);
    Free_ParamNames(all_vars, nvar+nparam);
    options_free(options);
}
