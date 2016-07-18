#include <stdio.h>
#include "fml.h"
#include "mag.h"
#include "rld0.h"
#include "htab.h"

int main(int argc, char *argv[])
{
	int n_seqs;
	bseq1_t *seqs;
	fml_opt_t opt;
	struct rld_t *e;
	struct mag_t *g;

	fml_opt_init(&opt);
	opt.bfc_opt.k = 17;
	if (argc == 1) return 1;
	seqs = bseq_read(argv[1], &n_seqs);
	fml_correct(&opt, n_seqs, seqs);
//	int i; for (i = 0; i < n_seqs; ++i) printf("@%s\n%s\n+\n%s\n", seqs[i].name, seqs[i].seq, seqs[i].qual);
	e = fml_fmi_gen(&opt, n_seqs, seqs);
//	rld_dump(e, "-");
	g = fml_assemble(&opt, e);
	fml_graph_clean(&opt, g);
	mag_g_print(g);
	return 0;
}
