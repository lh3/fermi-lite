#include <assert.h>
#include "internal.h"
#include "kstring.h"
#include "rle.h"
#include "mrope.h"
#include "rld0.h"
#include "mag.h"
#include "kvec.h"
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
	opt->min_merge_len = 0;
	bfc_opt_init(&opt->bfc_opt);
	mag_init_opt(&opt->mag_opt);
	opt->mag_opt.flag = MAG_F_CLEAN | MAG_F_NO_SIMPL;
	opt->bfc_opt.n_threads = opt->n_threads;
}

static inline int is_rev_same(int l, const char *s)
{
	int i;
	if (l&1) return 0;
	for (i = 0; i < l>>1; ++i)
		if (s[i] + s[l-1-i] != 5) break;
	return (i == l>>1);
}

struct rld_t *fml_fmi_gen(int n, bseq1_t *seq, int is_mt)
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
		free(s->qual); free(s->name);
		for (i = 0; i < s->l_seq; ++i)
			s->seq[i] = seq_nt6_table[(int)s->seq[i]];
		for (i = 0; i < s->l_seq; ++i)
			if (s->seq[i] == 5) break;
		if (i < s->l_seq) {
			free(s->seq);
			continue;
		}
		if (is_rev_same(s->l_seq, s->seq))
			--s->l_seq, s->seq[s->l_seq] = 0;
		seq_reverse(s->l_seq, (uint8_t*)s->seq);
		kputsn(s->seq, s->l_seq + 1, &str);
		seq_revcomp6(s->l_seq, (uint8_t*)s->seq);
		kputsn(s->seq, s->l_seq + 1, &str);
		free(s->seq);
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

struct rld_t *fml_seq2fmi(const fml_opt_t *opt, int n, bseq1_t *seq)
{
	return fml_fmi_gen(n, seq, opt->n_threads > 1? 1 : 0);
}

void fml_fmi_destroy(rld_t *e)
{
	rld_destroy(e);
}

void fml_mag_clean(const fml_opt_t *opt, struct mag_t *g)
{
	magopt_t o = opt->mag_opt;
	o.min_merge_len = opt->min_merge_len;
	mag_g_merge(g, 1, opt->min_merge_len);
	mag_g_clean(g, &o);
	mag_g_trim_open(g, &o);
}

void fml_mag_destroy(struct mag_t *g)
{
	mag_g_destroy(g);
}

#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

#define edge_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0) // from mag.c

fml_utg_t *fml_mag2utg(struct mag_t *g, int *n)
{
	size_t i, j;
	fml_utg_t *utg;
	khash_t(64) *h;
	khint_t k;

	h = kh_init(64);
	for (i = j = 0; i < g->v.n; ++i) {
		int absent;
		magv_t *p = &g->v.a[i];
		if (p->len < 0) continue;
		k = kh_put(64, h, p->k[0], &absent);
		kh_val(h, k) = j<<1 | 0;
		k = kh_put(64, h, p->k[1], &absent);
		kh_val(h, k) = j<<1 | 1;
		++j;
	}
	*n = j;
	kh_destroy(64, g->h);

	utg = (fml_utg_t*)calloc(*n, sizeof(fml_utg_t));
	for (i = j = 0; i < g->v.n; ++i) {
		magv_t *p = &g->v.a[i];
		fml_utg_t *q;
		int from, a, b;
		if (p->len < 0) continue;
		q = &utg[j++];
		q->len = p->len, q->nsr = p->nsr;
		q->seq = p->seq, q->cov = p->cov;
		for (from = 0; from < 2; ++from) {
			ku128_v *r = &p->nei[from];
			for (b = q->n_ovlp[from] = 0; b < r->n; ++b)
				if (!edge_is_del(r->a[b])) ++q->n_ovlp[from];
		}
		q->ovlp = (fml_ovlp_t*)calloc(q->n_ovlp[0] + q->n_ovlp[1], sizeof(fml_ovlp_t));
		for (from = a = 0; from < 2; ++from) {
			ku128_v *r = &p->nei[from];
			for (b = 0; b < r->n; ++b) {
				ku128_t *s = &r->a[b];
				fml_ovlp_t *t;
				if (edge_is_del(*s)) continue;
				t = &q->ovlp[a++];
				k = kh_get(64, h, s->x);
				assert(k != kh_end(h));
				t->tid = kh_val(h, k);
				t->len = s->y;
				t->from = from;
			}
			free(p->nei[from].a);
		}
	}
	kh_destroy(64, h);
	free(g->v.a);
	free(g);
	return utg;
}

void fml_utg_print(int n, const fml_utg_t *utg)
{
	int i, j, l;
	kstring_t out = {0,0,0};
	for (i = 0; i < n; ++i) {
		const fml_utg_t *u = &utg[i];
		out.l = 0;
		kputc('@', &out); kputw(i<<1|0, &out); kputc(':', &out); kputw(i<<1|1, &out);
		kputc('\t', &out); kputw(u->nsr, &out);
		kputc('\t', &out);
		for (j = 0; j < u->n_ovlp[0]; ++j) {
			kputw(u->ovlp[j].tid, &out); kputc(',', &out);
			kputw(u->ovlp[j].len, &out); kputc(';', &out);
		}
		if (j == 0) kputc('.', &out);
		kputc('\t', &out);
		for (; j < u->n_ovlp[0] + u->n_ovlp[1]; ++j) {
			kputw(u->ovlp[j].tid, &out); kputc(',', &out);
			kputw(u->ovlp[j].len, &out); kputc(';', &out);
		}
		kputc('\n', &out);
		l = out.l;
		kputsn(u->seq, u->len, &out);
		for (j = l; j < out.l; ++j)
			out.s[j] = "ACGTN"[(int)out.s[j] - 1];
		kputsn("\n+\n", 3, &out);
		kputsn(u->cov, u->len, &out);
		kputc('\n', &out);
		fputs(out.s, stdout);
	}
	free(out.s);
}

void fml_utg_destroy(int n, fml_utg_t *utg)
{
	int i;
	for (i = 0; i < n; ++i) {
		free(utg[i].seq);
		free(utg[i].cov);
		free(utg[i].ovlp);
	}
}

fml_utg_t *fml_assemble(const fml_opt_t *opt, int n_seqs, bseq1_t *seqs, int *n_utg)
{
	rld_t *e;
	mag_t *g;
	fml_utg_t *utg;

	fml_correct(opt, n_seqs, seqs);
	e = fml_seq2fmi(opt, n_seqs, seqs);
	g = fml_fmi2mag(opt, e);
	fml_mag_clean(opt, g);
	utg = fml_mag2utg(g, n_utg);
	return utg;
}
