#################
# Finding Genes #
#################

# Finding genes of interest in the genome assembly.

## Common stuff needed for all genes of interes

### Genome assembly sequence

SoybeanLooperScaffolds.fasta: ../genome/SoybeanLooperScaffolds.fasta
	ln -f ../genome/SoybeanLooperScaffolds.fasta ./

### Blast database of genome assembly

#### Note we use the database header file as the target

SoybeanLooperScaffolds.fasta.nhr: SoybeanLooperScaffolds.fasta
	makeblastdb -dbtype nucl -in SoybeanLooperScaffolds.fasta

.PHONY: blastdb

blastdb: SoybeanLooperScaffolds.fasta.nhr

## Ryanodine receptor

### tblastn search to find scaffolds containg putatuve ryanodine receptor exons

ryanodineReceptorHits.txt: SoybeanLooperScaffolds.fasta.nhr
	tblastn -db SoybeanLooperScaffolds.fasta -query noctuidRyanodineReceptor.fasta -evalue 1e-09 -outfmt 7 > ryanodineReceptorHits.txt

.PHONY: blastryanodine

blastryanodine: ryanodineReceptorHits.txt

### extract the putative exon-containing scaffolds to a new fasta file

ryanodineReceptorScaffs.txt: ryanodineReceptorHits.txt
	grep -v "#" ryanodineReceptorHits.txt | cut -f2 | sort -u > ryanodineReceptorScaffs.txt

ryanodineReceptorScaffs.fasta: ryanodineReceptorScaffs.txt
	seqkit grep -n -f ryanodineReceptorScaffs.txt SoybeanLooperScaffolds.fasta > ryanodineReceptorScaffs.fasta

.PHONY: ryanodinescaffolds

ryanodinescaffolds: ryanodineReceptorScaffs.fasta

### run exonerate with H. armigera as query

#### first, just get the "vulgar" lines and sort to get coverage

H.armigera.hits.sorted.txt: ryanodineReceptorScaffs.fasta H.armigera.ryanodineReceptor.fasta
	exonerate --model protein2genome --showquerygff yes \
	--showtargetgff yes H.armigera.ryanodineReceptor.fasta ryanodineReceptorScaffs.fasta \
	| grep "^vulgar" | sort -g -k3 | sed 's/vulgar: //' > H.armigera.hits.sorted.txt

#### comprehensive exonerate output including alignments and gtf blocks

H.armigera.ryanodine.exonerate.out: ryanodineReceptorScaffs.fasta H.armigera.ryanodineReceptor.fasta
	exonerate --model protein2genome --showquerygff yes \
	--showtargetgff yes H.armigera.ryanodineReceptor.fasta ryanodineReceptorScaffs.fasta \
	> H.armigera.ryanodine.exonerate.out

.PHONY: H.armigera_exonerate

H.armigera_exonerate: H.armigera.hits.sorted.txt H.armigera.ryanodine.exonerate.out

#### Pulling exon 1 and exon 1 contigs

H.armigera.exon1.fasta: exon1.gtf ryanodineReceptorScaffs.fasta
	seqkit subseq --gtf exon1.gtf ryanodineReceptorScaffs.fasta > exon1.fasta

H.armigera.exon1.scaffs.fasta: ryanodineReceptorScaffs.fasta
	seqkit grep -p scaffold42287,scaffold27611 ryanodineReceptorScaffs.fasta \
	> H.armigera.exon1.scaffs.fasta

.PHONY: exon1

exon1: H.armigera.exon1.fasta H.armigera.exon1.scaffs.fasta
