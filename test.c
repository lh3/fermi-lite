#include <stdio.h>
#include "fml.h"
#include "htab.h"

int main(int argc, char *argv[])
{
	int n_seqs;
	bseq1_t *seqs;
	fml_opt_t opt;

	fml_opt_init(&opt);
	if (argc == 1) return 1;
	seqs = bseq_read(argv[1], &n_seqs);
	fml_correct(&opt, n_seqs, seqs);
	return 0;
}
