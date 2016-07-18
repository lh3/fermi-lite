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
	float min_dratio1, max_bcov, max_bfrac;
} magopt_t;

typedef struct {
	int n_threads, min_match, min_merge_len;
	bfc_opt_t bfc_opt;
	magopt_t mag_opt;
} fml_opt_t;

struct rld_t;
struct mag_t;

typedef struct {
	uint32_t tid;
	uint32_t len:31, from:1;
} fml_ovlp_t;

typedef struct {
	int32_t len;      // length of sequence
	int32_t nsr;      // number of supporting reads
	char *seq;        // unitig sequence, in the nt6 encoding
	char *cov;        // per-base coverage
	int n_ovlp[2];    // number of 5'-end [0] and 3'-end [1] overlaps
	fml_ovlp_t *ovlp; // overlaps, of size n_ovlp[0]+n_ovlp[1]
} fml_utg_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize default parameters
 *
 * @param opt      (out) pointer to parameters
 */
void fml_opt_init(fml_opt_t *opt);

void fml_opt_set_ec_k(fml_opt_t *opt, int k);
void fml_opt_set_min_ovlp(fml_opt_t *opt, int min_ovlp);
void fml_opt_set_clean_ovlp(fml_opt_t *opt, int min_ovlp, float min_drop_ratio);
void fml_opt_set_merge_ovlp(fml_opt_t *opt, int min_merge_ovlp);

/**
 * Assemble a list of sequences
 *
 * @param opt      parameters
 * @param n_seqs   number of input sequences
 * @param seqs     sequences to assemble; FREED on return
 * @param n_utg    (out) number of unitigs in return
 *
 * @return list of unitigs
 */
fml_utg_t *fml_assemble(const fml_opt_t *opt, int n_seqs, bseq1_t *seqs, int *n_utg);

/**
 * Free unitigs
 *
 * @param n_utg    number of unitigs
 * @param utg      list of unitigs
 */
void fml_utg_destroy(int n_utg, fml_utg_t *utg);

/**
 * Read all sequences from FASTA/FASTQ
 *
 * @param fn       filename; NULL or "-" for stdin
 * @param n        (out) number of sequences read into RAM
 *
 * @return list of sequences
 */
bseq1_t *bseq_read(const char *fn, int *n);

/**
 * Error correction
 *
 * @param opt       parameters
 * @param n         number of sequences
 * @param seq       sequences; corrected IN PLACE
 */
void fml_correct(const fml_opt_t *opt, int n, bseq1_t *seq);

/**
 * Construct FMD-index
 *
 * @param opt       parameters
 * @param n         number of sequences
 * @param seq       sequences; FREED on return
 *
 * @return FMD-index
 */
struct rld_t *fml_seq2fmi(const fml_opt_t *opt, int n, bseq1_t *seq);

/**
 * Generate initial overlap graph
 *
 * @param opt       parameters
 * @param e         FMD-index; FREED on return
 *
 * @return overlap graph in the "mag" structure
 */
struct mag_t *fml_fmi2mag(const fml_opt_t *opt, struct rld_t *e);

/**
 * Clean a mag graph
 *
 * @param opt       parameters
 * @param g         overlap graph; modified IN PLACE
 */
void fml_mag_clean(const fml_opt_t *opt, struct mag_t *g);

/**
 * Convert a graph in mag to fml_utg_t
 *
 * @param g         graph in the "mag" structure; FREED on return
 * @param n_utg     (out) number of unitigs
 *
 * @return list of unitigs
 */
fml_utg_t *fml_mag2utg(struct mag_t *g, int *n_utg);

void fml_utg_print(int n, const fml_utg_t *utg);
void fml_fmi_destroy(struct rld_t *e);
void fml_mag_destroy(struct mag_t *g);

#ifdef __cplusplus
}
#endif

#endif
