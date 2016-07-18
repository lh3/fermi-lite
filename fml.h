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

#define MAG_F_READ_ORI   0x1
#define MAG_F_READ_TAG   0x2
#define MAG_F_CLEAN      0x10
#define MAG_F_AGGRESSIVE 0x20
#define MAG_F_NO_AMEND   0x40
#define MAG_F_NO_SIMPL   0x80

typedef struct {
	int flag, max_arc, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bvtx, min_merge_len, trim_len, trim_depth;
	float min_dratio0, min_dratio1;
	float max_bcov, max_bfrac;
} magopt_t;

typedef struct {
	int n_threads, min_match, min_merge_len;
	bfc_opt_t bfc_opt;
	magopt_t mag_opt;
} fml_opt_t;

struct rld_t;
struct mag_t;

#ifdef __cplusplus
extern "C" {
#endif

void fml_opt_init(fml_opt_t *opt);
bseq1_t *bseq_read(const char *fn, int *n_);
void fml_correct(const fml_opt_t *opt, int n, bseq1_t *seq);
struct rld_t *fml_fmi_gen(const fml_opt_t *opt, int n, bseq1_t *seq);
void fml_fmi_destroy(struct rld_t *e);
struct mag_t *fml_assemble(const fml_opt_t *opt, const struct rld_t *e);
void fml_graph_clean(const fml_opt_t *opt, struct mag_t *g);
void fml_graph_destroy(struct mag_t *g);

#ifdef __cplusplus
}
#endif

#endif
