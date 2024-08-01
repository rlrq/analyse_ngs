# analyse_ngs
Process NGS reads from paired fastq files & booking sheet, for chaelab.

## Steps
1. Dereplicate reads using dada2 (R)
2. BLASTN of unique read sequences to reference genome (FASTA)
3. Intersect BLAST hits with reference genome annotations (GFF) using bedtools
4. Map unique read sequences to best hit's gene ID ('.' for hits not overlapping with genes)
5. Summarise abundance of reads mapped to each gene per run
6. Demultiplex by separating reads by mapped gene using seqkit, retained only if genes are on-target per booking sheet
