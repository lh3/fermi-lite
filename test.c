#include <stdio.h>
#include "fml.h"
#include "htab.h"

int main(int argc, char *argv[])
{
	int i, n_seqs;
	bseq1_t *seqs;
	fml_opt_t opt;

	fml_opt_init(&opt);
	opt.bfc_opt.k = 17;
	if (argc == 1) return 1;
	seqs = bseq_read(argv[1], &n_seqs);
	fml_correct(&opt, n_seqs, seqs);
	for (i = 0; i < n_seqs; ++i)
		printf("@%s\n%s\n+\n%s\n", seqs[i].name, seqs[i].seq, seqs[i].qual);
	return 0;
}
