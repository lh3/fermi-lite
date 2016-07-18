#include "internal.h"
#include "fml.h"

void fml_opt_init(fml_opt_t *opt)
{
	bfc_opt_init(&opt->bfc_opt);
}
