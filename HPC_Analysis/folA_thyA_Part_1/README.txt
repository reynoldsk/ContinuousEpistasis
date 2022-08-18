README

Written 11/1/21 by Ryan Otto

This code was used to count sgRNA frequencies over time from pairwise CRISPRi HiSeq data. As noted below, the FASTQ files themselves are not present, but can be downloaded separately (see publication for details).

Contents:
BarcodedCounts -- folder containing pairwise count data for all timepoints
nt_sgRNA.txt -- Key indicating the critical nucleotides for determining sgRNA identity.
runPython -- Bash command to run the python counting script
sgRNAs_seq.fa -- fasta file containing nucleotide sequences for all sgRNAs
step2_code.py -- python file for deciphering HiSeq counts. Name of file is an artifact of previous analysis that involved a read merging step.
*_files.txt -- All of these are simple text files for parallelizing the counting setup.

Lacks:
store_fastq -- A folder that, for the code to run, must contain un-gzipped HiSeq FASTQ files.