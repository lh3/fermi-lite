#include "internal.h"
#include "kstring.h"
#include "rle.h"
#include "mrope.h"
#include "rld0.h"
#include "mag.h"
#include "fml.h"

unsigned char seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void fml_opt_init(fml_opt_t *opt)
{
	opt->n_threads = 1;
	opt->min_match = 55;
	opt->min_merge_len = 71;
	bfc_opt_init(&opt->bfc_opt);
	mag_init_opt(&opt->mag_opt);
}

static inline int is_rev_same(int l, const char *s)
{
	int i;
	if (l&1) return 0;
	for (i = 0; i < l>>1; ++i)
		if (s[i] + s[l-1-i] != 5) break;
	return (i == l>>1);
}

struct rld_t *fml_fmi_gen_core(int n, bseq1_t *seq, int is_mt)
{
	mrope_t *mr;
	kstring_t str = {0,0,0};
	mritr_t itr;
	rlditr_t di;
	const uint8_t *block;
	rld_t *e = 0;
	int k;

	mr = mr_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN, MR_SO_RCLO);
	for (k = 0; k < n; ++k) {
		int i;
		bseq1_t *s = &seq[k];
		for (i = 0; i < s->l_seq; ++i)
			s->seq[i] = seq_nt6_table[(int)s->seq[i]];
		for (i = 0; i < s->l_seq; ++i)
			if (s->seq[i] == 5) break;
		if (i < s->l_seq) continue;
		if (is_rev_same(s->l_seq, s->seq))
			--s->l_seq, s->seq[s->l_seq] = 0;
		kputsn(s->seq, s->l_seq + 1, &str);
		seq_revcomp6(s->l_seq, (uint8_t*)s->seq);
		if (s->l_seq&1) s->seq[i] = 5 - s->seq[i];
		kputsn(s->seq, s->l_seq + 1, &str);
		free(s->seq); free(s->qual); free(s->name);
	}
	free(seq);
	mr_insert_multi(mr, str.l, (uint8_t*)str.s, is_mt);
	free(str.s);

	e = rld_init(6, 3);
	rld_itr_init(e, &di, 0);
	mr_itr_first(mr, &itr, 1);
	while ((block = mr_itr_next_block(&itr)) != 0) {
		const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
		while (q < end) {
			int c = 0;
			int64_t l;
			rle_dec1(q, c, l);
			rld_enc(e, &di, l, c);
		}
	}
	rld_enc_finish(e, &di);

	mr_destroy(mr);
	return e;
}

struct rld_t *fml_fmi_gen(const fml_opt_t *opt, int n, bseq1_t *seq)
{
	return fml_fmi_gen_core(n, seq, opt->n_threads > 1? 1 : 0);
}

void fml_fmi_destroy(rld_t *e)
{
	rld_destroy(e);
}

void fml_graph_clean(const fml_opt_t *opt, struct mag_t *g)
{
	magopt_t o = opt->mag_opt;
	o.min_merge_len = opt->min_merge_len;
	mag_g_merge(g, 1, opt->min_merge_len);
	mag_g_clean(g, &o);
	mag_g_trim_open(g, &o);
}

void fml_graph_destroy(struct mag_t *g)
{
	mag_g_destroy(g);
}
