# Soybean looper MIPS analysis

Note that prior to starting analysis, I set up a fresh conda environment (SoybeanLooper) to install needed programs via bioconda wherever possible. Analysis should be run from this environment. Steps in the analysis are recorded, and can be repreated using the Makefil in this directory. The Makefile, these notes and other small text files are kept under a git repository. Large files are not kept in the repository. In particular, fastq files and other output produced by the MiniSeq are in a directory <MiniSeqOutput> that is not tracked by git.

## Generating genomic sequence

To generate some genomic contigs, sequenced a single individual: MA22 with one lane of MiniSeq. Original data are stored in <MiniSeqOutPut/20180331_155401>. To keep things easy to follow, set up a directory <GenomeRawData> and linked the fastq files into it.

### Initial data QC.

Ran fastqc 0.11.7 on the two raw data files. Total number of read pairs is 28,001,390. Overall, the data look decent. Theere is as usual, both a decline in mean quality score and an increase in variance in quality score toward the end of the reads. Fastqc does not find any adaptor sequences. There is one overrepresented sequence in the read 2 data, but not iditifiable source, and it is found in 0.19% of reads.

Conclusion - we are in good shape, but we will want to do some quality-trimming before attempting assembly.
