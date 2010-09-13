#include <isl_obj.h>

struct isl_str {
	int ref;

	struct isl_ctx *ctx;

	char *s;
};
typedef struct isl_str		isl_str;

__isl_give isl_str *isl_str_alloc(struct isl_ctx *ctx);
__isl_give isl_str *isl_str_concat(__isl_take isl_str *str1,
	__isl_take isl_str *str2);

extern struct isl_obj_vtable isl_obj_str_vtable;
#define isl_obj_str		(&isl_obj_str_vtable)