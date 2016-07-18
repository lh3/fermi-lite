#ifndef FML_H
#define FML_H

#include <stdint.h>

typedef struct {
	int l_seq;
	char *name, *seq, *qual;
} bseq1_t;

typedef struct {
	int n_threads, q, k, l_pre;
	int min_cov; // a k-mer is considered solid if the count is no less than this

	int max_end_ext;
	int win_multi_ec;

	// these ec options cannot be changed on the command line
	int w_ec, w_ec_high, w_absent, w_absent_high;
	int max_path_diff, max_heap;
} bfc_opt_t;

typedef struct {
	bfc_opt_t bfc_opt;
} fml_opt_t;

#ifdef __cplusplus
extern "C" {
#endif

void fml_opt_init(fml_opt_t *opt);
bseq1_t *bseq_read(const char *fn, int *n_);
void fml_correct(const fml_opt_t *opt, int n, bseq1_t *seq);

#ifdef __cplusplus
}
#endif

#endif
