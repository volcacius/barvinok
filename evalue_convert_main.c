#include "argp.h"
#include <barvinok/options.h>
#include "evalue_convert.h"
#include "evalue_read.h"

int main(int argc, char **argv)
{
    evalue *EP;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    struct convert_options   convert;
    int printed;

    set_program_name(argv[0]);
    argp_parse(&convert_argp, argc, argv, 0, 0, &convert);

    EP = evalue_read_from_file(stdin, NULL, &all_vars,
			       &nvar, &nparam, options->MaxRays);
    assert(EP);

    printed =
	evalue_convert(EP, &convert, options->verbose, nvar+nparam, all_vars);

    if (!printed)
	print_evalue(stdout, EP, all_vars);

    evalue_free(EP);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(options);
}
