#include <omega.h>
#include <barvinok/util.h>
#include "omega/convert.h"
#include "vertices.h"
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

void vertices(Relation& r)
{
    varvector vv;
    varvector params;
    Param_Polyhedron *PP;

    Polyhedron *D = relation2Domain(r, vv, params);
    assert(!D->next);
    Polyhedron *ctx = Universe_Polyhedron(params.size());

    r.setup_names();
    const char **param_names = new const char *[params.size()];
    for (int i = 0; i < params.size(); ++i)
	param_names[i] = params[i]->char_name();

    PP = Polyhedron2Param_Domain(D, ctx, MAXRAYS);
    Param_Polyhedron_Print(stdout, PP, (char **)param_names);

    delete [] param_names;

    Param_Polyhedron_Free(PP);
    Polyhedron_Free(ctx);
    Domain_Free(D);
}
