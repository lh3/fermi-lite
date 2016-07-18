#ifndef FML_INTERNAL_H
#define FML_INTERNAL_H

#include "fml.h"

extern unsigned char seq_nt6_table[256];

#ifdef __cplusplus
extern "C" {
#endif

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void bfc_opt_init(bfc_opt_t *opt);
void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

#ifdef __cplusplus
}
#endif

#endif
