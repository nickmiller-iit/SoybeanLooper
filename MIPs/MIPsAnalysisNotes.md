## First test run

The first test sequencing run consisted of a total of 64 individauls. Ten were fall armworms, 46 were soybean loopers with 8 soybean loopers run in duplicate. Sequencing was succesful insofar as I got data and the indexing reads split the data into individual samples.

On a side note, I'm trying out snakemake as an alternative to make. It's got a bit of learning curve but it looks like it will work better running jobs in parallel when there are multiple input or output files per job.

### Initial QC

### Trimming

Reads are 151 bases long. Inserts range from 125 - 138 bp, so we definitely expect to see primer sequences showing up at the end of the reads. Used trimmomatic to remove the primer sequences. Because the primers inlude index sequences, set up the file with just the parts of the primers 3' to the indexes.

Confusingly, I don't seem to be finding any primer sequence.

Explanation - I was wrong about the insert size, because it include the extension and ligation arms. All insers are 154 bp, so we do not expect to see any primer.

### Read merging

Reads should almost span the entire insert. We can merge read pairs into singe leads with PEAR. Looking at the fastqc results for the trimmed reads there was a small population of reads with lengths <= 80 bases. Re-ran trimmomatic with minlen increased to 85 to get rid of these

### Alignment

Ran alignment with bwa mem and default parameters. Visualized aligned reads with Tablet. I see three basic patterns.

1. Loci where no reads align
2. Loci where reads align along their full lengths. Eyeballing suggests plenty of heterozygotes
3. Loci where reads only align along part of their length and are "soft clipped"

This third category is a bit of a worry. Some quick eye\balling indicates that the clipped reads match to the target locus at both ends but not in the middle. This indicates that they are off-target reads from genome locations that match the probe arms. This could be a result of repetitive elements, despite our best efforts to avoid those.

In principle, a simple way to deal with this is to not include / get rid of soft-clipped reads. Bwa pretty much insists on soft clipping but bowtie, for example defaults to end-to-end alignments. Using bowtie 2 does indeed get rid of the clipped reads. This results in some additional loci wuth no reads aligning. That's OK, because they were not useable loci in the first place.

A comparison with MIPs data from fall armyworm partially supports the hypothesis that the off-target clipped reads were due to repetitive elements. The FAW MIPs were designed from a published genome sequence with repeats masked. I see a lot fewer loci in FAW with soft-clipped reads when aligning tith bwa. Nevertheless there are still a couple of loci with clipped reads, so repeats may no be the only issue

### Genotype calling and filtering

Initially, I thought I'd try treating the MIPs like genome-mapped RAD-Seq reads (there are some similarities after all). Trying stacks was a bust because it is finnicky about almost everything. Looking at the docs, dDocent just uses freebayes to call polymorphisms, so we might as well use freebayes directly.

Freebayes calls a lot of garbage polymorphisms, so we then have to filter out the stuff that is likely to be bogus. The first step was to filter out non-SNP polymorphisms. Theres aren't neccessarily likely to be bogus, but they are tough to model. The next step is to drop polymorphisms with low overall quality scores or highly unbalanced "alleles" in the heterozygotes. Thie third filter gets rid of any indivdiuals for which a majority of the polymorphisms are not genotyped and any polymophisms for which less than 80% of the individuals are genotyped.

These initial filters leave us with 39 SNPs distributed over just 7 loci. For most purposes, having multiple SNPs withn a few 10s of bp of each other is not terrible helpful because they are likely to be in LD. So we have a couple of options

 * Pick one SNP per locus
 * Try to resolve the SNPs into haplotypes.
 
Attempting to resolve haplotypes with freebayes doesn't give us any haplotypes. 

Turns out that most of the SNPs are rare, only a handful at intermediate frequencies.
