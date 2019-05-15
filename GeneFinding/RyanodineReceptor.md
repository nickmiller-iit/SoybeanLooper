# Finding the ryanodine receptor

A test case - how much of a gene of interest can we recover from the assembly. Start with the ryanodine receptor, the target for diamides.

## Example ryo receptors

Downloaded 5 full length protein sequences for ryanodin receptors from other noctuids"

 * *Spodoptera exigua*
 * *Helicoverpa armigera*
 * *Sesamia inferens*
 * *Mythimna separata*
 * *Spodoptera litura*

A quick alignment in Genious shows the aa sequence is pretty highly conserved. Total length of the aa sequence ranges from 5098 to 5142.

## Searching the genome assembly with tblastn

Using protein sequences from other noctuids to search the nucleotide genome sequence should get us candidate contigs containing parts of the ryanodine receptor.

Searching with an e-value cutoff of 10^-9, we get 49 unique scaffolds with putative exons.

## Finding exons with exonerate

Now that we've fished out putative exon-bearing contigs, use exonerate to align protein sequences to genome. It gets a bit messy if we align multiple protein sequences, so just use *H. armigera*, the longest one. Should be OK, given the protein is pretty conserved.

First, take the "vulgar" lines from the exonerate output and sort by start position of protein hit. The idea here is to make it east to see how much of the protein we have covered and where there are gaps.  At this point, things got pretty hairy. There are a whole bunch of cases of the same part of the protein hitting to different scaffolds.

AHB33498.1 15 87 . scaffold27611 3516 1259
AHB33498.1 19 74 . scaffold42287 1693 1858

These two have perfect hits from aa 20 to 70. 27611 has a 2041bp intron, then a short exon up to aa 87

Extracting the first exon as "exon 1" and the 2 contigs as exon 1 contigs, we can align them and see what's up using Genious. The nucleotide sequences of the exonic regions are almost identical, 2 synonymous differences. The flanking regions, however, are very different.

----

AHB33498.1 92 428 . scaffold10712 1378 6668 +
AHB33498.1 133 252 . scaffold30469 3780 1121 -

scaffold10712 has a hit split into 7 exons. Let's call the first of these exon 3, since it more or less picks up from exon 2 of scaffold27611.

scaffold30469 has a hit split into 3 exons. These appear to correspond (roughly) to exons 4 (partial), 5, and 6 on scaffold10712.

Again, if we align everything, the exons are very similar, only synonymous differences. Th introns however are very different.

My suspicion at this stage is we are dealing with a duplicate gene. The alternative is that we have alleles assembling as spearate contigs. This is possible, but if it is common, we would expect a lot of duplicated BUSCOs, which we don't see.

## How duplicated are the exons

It would be hopelessly time consuming to fo through all of the hits individually. An alterantive solution is:

 1. Extract just the exonic sequences corresponding to hits.
 2. Cluster exonic sequences with cd-hit
 3. Look at the distibution of cluster sizes
 
If we are dealing with a duplicate gene, we would expect most exons to show up twice.

It turns out every cluster contains either 1 or 2 sequences. The format of the cluster file is not especiaaly easy to parse into counts of sequences per cluster but it can be done manually with grep, cu, sort and wc.

+----------+----------+
|Sequences |Count     |
+----------+----------+
|1         |107       |
+----------+----------+
|2         |26        |
+----------+----------+

So, maybe it's not a duplicated gene after all. Overall, we need better scaffolding to tell.

