#ifndef FML_H
#define FML_H

#include <stdint.h>

typedef struct {
	int l_seq;
	uint32_t aux, aux2;
	char *name, *comment, *seq, *qual;
} bseq1_t;

typedef struct {
	int n_threads, q, k, l_pre;
	float min_frac;

	int max_end_ext;
	int win_multi_ec; // no 2 high-qual corrections or 4 corrections in a window of this size
	int min_cov; // a k-mer is considered solid if the count is no less than this

	// these ec options cannot be changed on the command line
	int w_ec, w_ec_high, w_absent, w_absent_high;
	int max_path_diff, max_heap;
} bfc_opt_t;

typedef struct {
	bfc_opt_t bfc_opt;
} fml_opt_t;

struct bfc_ch_s;

#ifdef __cplusplus
extern "C" {
#endif

void fml_opt_init(fml_opt_t *opt);
bseq1_t *bseq_read(const char *fn, int *n_);
struct bfc_ch_s *fml_count(const fml_opt_t *opt, int n, bseq1_t *seq);

#ifdef __cplusplus
}
#endif

#endif
