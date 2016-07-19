#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "fml.h"

int main(int argc, char *argv[])
{
	fml_opt_t opt;
	int c, n_seqs, n_utg;
	bseq1_t *seqs;
	fml_utg_t *utg;

	fml_opt_init(&opt);
	while ((c = getopt(argc, argv, "Ae:l:r:t:c:")) >= 0) {
		if (c == 'e') opt.ec_k = atoi(optarg);
		else if (c == 'c') opt.ec_min_cov = atoi(optarg);
		else if (c == 'l') opt.min_asm_ovlp = atoi(optarg);
		else if (c == 'r') opt.mag_opt.min_dratio1 = atof(optarg);
		else if (c == 'A') opt.mag_opt.flag |= MAG_F_AGGRESSIVE;
		else if (c == 't') opt.n_threads = atoi(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fml-asm [options] <in.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -e INT     k-mer length for error correction (0 for auto; -1 to disable) [%d]\n", opt.ec_k);
		fprintf(stderr, "  -c INT     min occurrence for a k-mer to be considered solid in ec [%d]\n", opt.ec_min_cov);
		fprintf(stderr, "  -l INT     min overlap length [%d]\n", opt.min_asm_ovlp);
		fprintf(stderr, "  -r FLOAT   drop shorter overlaps [%g]\n", opt.mag_opt.min_dratio1);
		fprintf(stderr, "  -t INT     number of threads (don't use for small data sets) [%d]\n", opt.n_threads);
		fprintf(stderr, "  -A         discard heterozygotes\n");
		return 1;
	}
	seqs = bseq_read(argv[optind], &n_seqs);
	utg = fml_assemble(&opt, n_seqs, seqs, &n_utg);
	fml_utg_print(n_utg, utg);
	fml_utg_destroy(n_utg, utg);
	return 0;
}
