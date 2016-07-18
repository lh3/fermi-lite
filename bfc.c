#include <stdlib.h>
#include <string.h>
#include "htab.h"
#include "kmer.h"
#include "internal.h"
#include "fml.h"

/*******************
 *** BFC options ***
 *******************/

void bfc_opt_init(bfc_opt_t *opt)
{
	memset(opt, 0, sizeof(bfc_opt_t));
	opt->n_threads = 1;
	opt->q = 20;
	opt->k = 33;
	opt->l_pre = 16;

	opt->min_frac = .9;

	opt->min_cov = 3;
	opt->win_multi_ec = 10;
	opt->max_end_ext = 5;

	opt->w_ec = 1;
	opt->w_ec_high = 7;
	opt->w_absent = 3;
	opt->w_absent_high = 1;
	opt->max_path_diff = 15;
	opt->max_heap = 100;
}

/**********************
 *** K-mer counting ***
 **********************/

#define CNT_BUF_SIZE 256

typedef struct { // cache to reduce locking
	uint64_t y[2];
	int is_high;
} insbuf_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	const bfc_opt_t *opt;
	bfc_ch_t *ch;
	int *n_buf;
	insbuf_t **buf;
} cnt_step_t;

bfc_kmer_t bfc_kmer_null = {{0,0,0,0}};

static int bfc_kmer_bufclear(cnt_step_t *cs, int forced, int tid)
{
	int i, k, r;
	if (cs->ch == 0) return 0;
	for (i = k = 0; i < cs->n_buf[tid]; ++i) {
		r = bfc_ch_insert(cs->ch, cs->buf[tid][i].y, cs->buf[tid][i].is_high, forced);
		if (r < 0) cs->buf[tid][k++] = cs->buf[tid][i];
	}
	cs->n_buf[tid] = k;
	return k;
}

static void bfc_kmer_insert(cnt_step_t *cs, const bfc_kmer_t *x, int is_high, int tid)
{
	int k = cs->opt->k;
	uint64_t y[2], hash;
	hash = bfc_kmer_hash(k, x->x, y);
	if (bfc_ch_insert(cs->ch, y, is_high, 0) < 0) {
		insbuf_t *p;
		if (bfc_kmer_bufclear(cs, 0, tid) == CNT_BUF_SIZE)
			bfc_kmer_bufclear(cs, 1, tid);
		p = &cs->buf[tid][cs->n_buf[tid]++];
		p->y[0] = y[0], p->y[1] = y[1], p->is_high = is_high;
	}
}

static void worker_count(void *_data, long k, int tid)
{
	cnt_step_t *cs = (cnt_step_t*)_data;
	bseq1_t *s = &cs->seqs[k];
	const bfc_opt_t *o = cs->opt;
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	uint64_t qmer = 0, mask = (1ULL<<o->k) - 1;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(o->k, x.x, c);
			qmer = (qmer<<1 | (s->qual == 0 || s->qual[i] - 33 >= o->q)) & mask;
			if (++l >= o->k) bfc_kmer_insert(cs, &x, (qmer == mask), tid);
		} else l = 0, qmer = 0, x = bfc_kmer_null;
	}
}

struct bfc_ch_s *fml_count(const fml_opt_t *opt, int n, bseq1_t *seq)
{
	int i;
	cnt_step_t cs;
	cs.opt = &opt->bfc_opt, cs.n_seqs = n, cs.seqs = seq;
	cs.ch = bfc_ch_init(cs.opt->k, cs.opt->l_pre);
	cs.n_buf = calloc(cs.opt->n_threads, sizeof(int));
	cs.buf = calloc(cs.opt->n_threads, sizeof(void*));
	for (i = 0; i < cs.opt->n_threads; ++i)
		cs.buf[i] = malloc(CNT_BUF_SIZE * sizeof(insbuf_t));
	if (cs.opt->n_threads == 1) {
		for (i = 0; i < cs.n_seqs; ++i)
			worker_count(&cs, i, 0);
	} else kt_for(cs.opt->n_threads, worker_count, &cs, cs.n_seqs);
	for (i = 0; i < cs.opt->n_threads; ++i) free(cs.buf[i]);
	free(cs.buf); free(cs.n_buf);
	return cs.ch;
}
