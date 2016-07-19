## Introduction

Fermi-lite is a standalone C *library* for assembling Illumina short reads in
regions from 100bp to 10 million bp in size. Unlike mainstream *de novo*
assemblers which try to achieve good contiguity, fermi-lite aims to retain
heterozygous events. It can be used as a local reassembler for variant calling.

Fermi-lite integrates the C source code of [bfc][bfc] for error correction,
[ropebwt2][rb2] for FM-index construction and [fermi2][fm2] for overlap-based
assembly. It is largely an in-memory light-weight version of [fermikit][fk]
without generating any intermediate files on disk, and inherits the
performance and relatively small memory footprint of all these separate
projects.

## Usage

For now, see [example.c][example] for the basic use of the library. Here is a
sketch of the example:
```cpp
int i, n_seqs, n_utgs;
bseq1_t *seqs;
fml_utg_t *utgs;
seqs = bseq_read(fastq_file_name, &n_seqs); // or fill the array with callers' functions
fml_opt_init(&opt);
utgs = fml_assemble(&opt, n_seqs, seqs, &n_utgs);
for (i = 0; i < n_utgs; ++i)
	printf("@%d\n%s\n+\n%s\n", i+1, utgs[i].seq, utgs[i].cov);
fml_utg_destroy(n_utgs, utg);
```
The output is in fact a graph. You may have a look at the [header file][header]
for details.

## Limitations

1. Fermi-lite can assemble bacterial genomes. It is fairly fast on multiple
   threads. The continguity and accuracy is expected to be similar to the older
   assemblers, but may not compete with more recent mainstream assemblers.

2. Fermi-lite does not work with genomes more than tens of megabases. It would
   take too much memory as it stages all data in memory. For large genomes,
   please use [fermikit][fk] instead.

3. This is the first iteration of fermi-lite. It is still immarture. In
   particular, I hope fermi-lite can be smarter to automatically set various
   parameters based on input, which is very challenging given the high
   variability of input data.

[bfc]: http://github.com/lh3/bfc
[rb2]: http://github.com/lh3/ropebwt2
[fm2]: http://github.com/lh3/fermi2
[fk]: http://github.com/lh3/fermikit
[example]: https://github.com/lh3/fermi-lite/blob/master/example.c
[header]: https://github.com/lh3/fermi-lite/blob/master/fml.h
