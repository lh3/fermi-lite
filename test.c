#include <stdio.h>
#include "fml.h"
#include "htab.h"

int main(int argc, char *argv[])
{
	int n_seqs;
	bseq1_t *seqs;
	fml_opt_t opt;
	struct bfc_ch_s *h;

	fml_opt_init(&opt);
	if (argc == 1) return 1;
	seqs = bseq_read(argv[1], &n_seqs);
	h = fml_count(&opt, n_seqs, seqs);
	printf("%llu\n", bfc_ch_count(h));
	return 0;
}
