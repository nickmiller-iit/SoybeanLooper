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
