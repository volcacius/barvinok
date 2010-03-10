#include <isl_obj.h>

struct isl_list {
	int ref;

	struct isl_ctx *ctx;

	int n;
	struct isl_obj obj[1];
};

struct isl_list *isl_list_alloc(struct isl_ctx *ctx, int n);
void isl_list_free(struct isl_list *list);

extern struct isl_obj_vtable isl_obj_list_vtable;
#define isl_obj_list		(&isl_obj_list_vtable)
