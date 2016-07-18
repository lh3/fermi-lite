#ifndef BFC_BSEQ_H
#define BFC_BSEQ_H

#include <stdint.h>

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

typedef struct {
	int l_seq;
	uint32_t aux, aux2;
	char *name, *comment, *seq, *qual;
} bseq1_t;

#ifdef __cplusplus
extern "C" {
#endif

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
bseq1_t *bseq_read(bseq_file_t *fp, int chunk_size, int keep_comment, int *n_);

#ifdef __cplusplus
}
#endif

#endif
