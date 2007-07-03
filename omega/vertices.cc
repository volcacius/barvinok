#include <omega.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "omega/convert.h"
#include "vertices.h"
#include "param_util.h"
#include "config.h"

void vertices(Relation& r)
{
    varvector vv;
    varvector params;
    Param_Polyhedron *PP;
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    Polyhedron *D = relation2Domain(r, vv, params);
    assert(!D->next);
    Polyhedron *ctx = Universe_Polyhedron(params.size());

    r.setup_names();
    const char **param_names = new const char *[params.size()];
    for (int i = 0; i < params.size(); ++i)
	param_names[i] = params[i]->char_name();

    PP = Polyhedron2Param_Polyhedron(D, ctx, options);
    Param_Polyhedron_Print(stdout, PP, (char **)param_names);

    delete [] param_names;

    Param_Polyhedron_Free(PP);
    Polyhedron_Free(ctx);
    Domain_Free(D);
    barvinok_options_free(options);
}
