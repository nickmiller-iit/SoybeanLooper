#############################################################
#############################################################
## Master Makefile to map reads to loci and call genotypes ##
#############################################################
#############################################################


###############
# Basic setup #
###############

RAW_READ_DIR=raw_fastq

RAW_R1=$(shell ls $(RAW_READ_DIR)/*_R1_*.fastq.gz)
RAW_R2=$(shell ls $(RAW_READ_DIR)/*_R2_*.fastq.gz)
RAW_ALL=$(RAW_R1) $(RAW_R2)

# We can extract sample names from the read 1 file names provided the following hold true
# 1. Filenames begin with the sample name
# 2. The sample name is separated from the rest of the file name by an underscore
# 3. Sample names do not contain underscores

SAMPLES=$(shell echo $(notdir $(RAW_R1)) | tr " " "\n" | cut -d "_" -f1)

# A custom function to get the sample name from a file name
define getsample
	$(basename $(strip $(shell echo $(notdir $(1)) | cut -d "_" -f1)))
endef

###################
# QC and trimming #
###################

# QC of raw data

FASTQC_DIR=fastqc_out

$(FASTQC_DIR):
	if [ ! -d $(FASTQC_DIR) ]; then mkdir $(FASTQC_DIR); fi

RAW_FASTQC_DIR=$(FASTQC_DIR)/raw

$(RAW_FASTQC_DIR): $(FASTQC_DIR)
	if [ ! -d $(RAW_FASTQC_DIR) ]; then mkdir $(RAW_FASTQC_DIR); fi

RAW_FASTQC_OUTS=$(subst .fastq.gz,_fastqc.html,$(subst $(RAW_READ_DIR), $(RAW_FASTQC_DIR), $(RAW_ALL)))

$(RAW_FASTQC_OUTS): $(RAW_ALL) | $(RAW_FASTQC_DIR)
	fastqc -o $(RAW_FASTQC_DIR) $?

RAW_FASTQC_SUMMARY=$(RAW_FASTQC_DIR)/summary.html

$(RAW_FASTQC_SUMMARY): $(RAW_FASTQC_OUTS)
	multiqc -n $(RAW_FASTQC_SUMMARY) $(RAW_FASTQC_DIR)

.PHONY: qc_raw

qc_raw: $(RAW_FASTQC_SUMMARY)


# Trimming
# This gets a bit tricky. The preferred option is to use automatic variables. That should let us
# run jobs in parallel and update only what is needed as we add new data

TRIM_DIR=trimmed_fastq

$(TRIM_DIR):
	if [ ! -d $(TRIM_DIR) ]; then mkdir $(TRIM_DIR); fi

1P_TRIM=$(addprefix $(TRIM_DIR)/, $(addsuffix _1P.fastq.gz, $(SAMPLES)))
2P_TRIM=$(addprefix $(TRIM_DIR)/, $(addsuffix _2P.fastq.gz, $(SAMPLES)))
1UP_TRIM=$(addprefix $(TRIM_DIR)/, $(addsuffix _1UP.fastq.gz, $(SAMPLES)))
2UP_TRIM=$(addprefix $(TRIM_DIR)/, $(addsuffix _2UP.fastq.gz, $(SAMPLES)))

TRIM_ALL=$(1P_TRIM) $(2P_TRIM) $(1UP_TRIM) $(2UP_TRIM)

TRIM_OPTS=ILLUMINACLIP:primers.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:85

# For reasons that I don't understand using the filter function to retrieve
# a filename from a set of files fails if the directory name is included.
$(1P_TRIM): $(RAW_R1) | $(TRIM_DIR)
	trimmomatic PE \
	$(RAW_READ_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(RAW_R1))) \
	$(RAW_READ_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(RAW_R2))) \
	$(TRIM_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(1P_TRIM))) \
	$(TRIM_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(1UP_TRIM))) \
	$(TRIM_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(2P_TRIM))) \
	$(TRIM_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(2UP_TRIM))) \
	$(TRIM_OPTS)

.PHONY: trim

trim: $(1P_TRIM)

# Post-trimming QC

TRIM_FASTQC_DIR=$(FASTQC_DIR)/trimmed

$(TRIM_FASTQC_DIR): | $(FASTQC_DIR)
	if [ ! -d $(TRIM_FASTQC_DIR) ]; then mkdir $(TRIM_FASTQC_DIR); fi

TRIM_FASTQC_OUTS=$(subst .fastq.gz,_fastqc.html,$(subst $(TRIM_DIR), $(TRIM_FASTQC_DIR), $(TRIM_ALL)))

$(TRIM_FASTQC_OUTS): $(TRIM_ALL) | $(TRIM_FASTQC_DIR)
	fastqc -o $(TRIM_FASTQC_DIR) $?

TRIM_FASTQC_SUMMARY=$(TRIM_FASTQC_DIR)/summary.html

$(TRIM_FASTQC_SUMMARY): $(TRIM_FASTQC_OUTS)
	multiqc -n $(TRIM_FASTQC_SUMMARY) $(TRIM_FASTQC_DIR)

.PHONY: qc_trim

qc_trim: $(TRIM_FASTQC_SUMMARY)


######################
# Merge paired reads #
######################

# Combine paired reads into merged reads with PEAR

MERGED_DIR=merged_fastq

$(MERGED_DIR):
	if [ ! -d $(MERGED_DIR) ]; then mkdir $(MERGED_DIR); fi

MERGED_READS=$(addprefix $(MERGED_DIR)/, $(addsuffix _merged.fastq.gz,$(SAMPLES)))

#PEAR likes to tack on a bunch of extra extrensions and does not compress output
$(MERGED_READS):$(TRIM_ALL) | $(MERGED_DIR)
	pear \
	-f $(TRIM_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(1P_TRIM))) \
	-r $(TRIM_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(2P_TRIM))) \
	-o $@;
	gzip -c $(addsuffix .assembled.fastq,$@) > $@;
	rm $(addsuffix .*.fastq,$@);

.PHONY: merge

merge: $(MERGED_READS)

# QC Merged data

MERGED_FASTQC_DIR=$(FASTQC_DIR)/merged

$(MERGED_FASTQC_DIR): | $(FASTQC_DIR)
	if [ ! -d $(MERGED_FASTQC_DIR) ]; then mkdir $(MERGED_FASTQC_DIR); fi


MERGED_FASTQC_OUTS=$(subst .fastq.gz,_fastqc.html,$(subst $(MERGED_DIR), $(MERGED_FASTQC_DIR), $(MERGED_READS)))

$(MERGED_FASTQC_OUTS): $(MERGED_READS) | $(MERGED_FASTQC_DIR)
	fastqc -o $(MERGED_FASTQC_DIR) $?

MERGED_FASTQC_SUMMARY=$(MERGED_FASTQC_DIR)/summary.html

$(MERGED_FASTQC_SUMMARY): $(MERGED_FASTQC_OUTS)
	multiqc -n $(MERGED_FASTQC_SUMMARY) $(MERGED_FASTQC_DIR)

.PHONY: qc_merged

qc_merge: $(MERGED_FASTQC_SUMMARY)



###################################
# Generating a reference for loci #
###################################

# We don't need to map reads to the whole genome, only the loci captured by the
# MIPS. Conveniently, the first part of the first part of the "X.mip_key" field
# of the mipgen output (up to the first "/") is in the same format as samtools
# uses to specify genome regions.

# Get the genome scaffold and selected MIPgen output from the genome dir.

GENOME_FASTA=SoybeanLooperScaffolds.fasta

$(GENOME_FASTA): ../genome/$(GENOME_FASTA)
	ln -s ../genome/$(GENOME_FASTA) .

CHOSEN_MIPS=SoybeanLooper.Selected.MIPs.txt

$(CHOSEN_MIPS): ../genome/$(CHOSEN_MIPS)
	ln -s ../genome/$(CHOSEN_MIPS)  $(CHOSEN_MIPS)

#faidx of genome file
GENOME_FASTA_IDX=SoybeanLooperScaffolds.fasta.fai

$(GENOME_FASTA_IDX): $(GENOME_FASTA)
	samtools faidx $(GENOME_FASTA)

# create an alignments dir with fasta of just the target loci

ALIGN_DIR=aligned

$(ALIGN_DIR):
	if [ ! -d $(ALIGN_DIR) ]; then mkdir $(ALIGN_DIR); fi

MIPS_FASTA=$(ALIGN_DIR)/loci.fasta

$(MIPS_FASTA): $(GENOME_FASTA_IDX) $(CHOSEN_MIPS) | $(ALIGN_DIR)
	samtools faidx $(GENOME_FASTA) `tail -n50 $(CHOSEN_MIPS) | cut -d "/" -f1` > $(MIPS_FASTA)

MIPS_FASTA_IDX=$(MIPS_FASTA).fai

$(MIPS_FASTA_IDX): $(MIPS_FASTA)
	samtools faidx $(MIPS_FASTA)

.PHONY: mips_loci_ref

mips_loci_ref: $(MIPS_FASTA_IDX)

#####################################
# Align sample reads to target loci #
#####################################

# Set up bowtie index
# There are multiple index files, any can serve as the target

BOWTIE_IDX=$(MIPS_FASTA).1.bt2

$(BOWTIE_IDX): $(MIPS_FASTA)
	bowtie2-build $(MIPS_FASTA) $(MIPS_FASTA)

#align sample reads to reference

ALIGNMENTS=$(addprefix $(ALIGN_DIR)/,$(addsuffix .bam, $(SAMPLES)))

$(ALIGNMENTS): $(BOWTIE_IDX) $(MERGED_READS)
	bowtie2 \
	-x $(MIPS_FASTA) \
	-U $(MERGED_DIR)/$(filter $(call getsample,$@)_%, $(notdir $(MERGED_READS))) \
	--rg-id $(strip $(call getsample,$@)) \
	--rg SM:$(strip $(call getsample,$@)) \
	| samtools view -b > $@

SORTED_ALIGNMENTS=$(subst .bam,.sorted.bam,$(ALIGNMENTS))

$(SORTED_ALIGNMENTS): $(ALIGNMENTS)
	samtools sort -o $@ \
	$(subst .sorted.bam,.bam,$@)

SORTED_ALIGNMENTS_IDX=$(addsuffix .bai,$(SORTED_ALIGNMENTS))

$(SORTED_ALIGNMENTS_IDX): $(SORTED_ALIGNMENTS)
	samtools index $(subst .bai,,$@)

.PHONY: align

align: $(SORTED_ALIGNMENTS_IDX)

############################################
# Calling and filtering polymorphisms      #
#                                          #
# NB. We confine ourselves to SNPs, Indels #
# etc. are ignored                         #
############################################

VCF_DIR=vcf

$(VCF_DIR):
	if [ ! -d $(VCF_DIR) ]; then mkdir $(VCF_DIR); fi


# We only want to call SNPs in the "scanned" region, not the ends corrsponding to the probe arms
# freebayes lets us specify targets in BED format. The required  fields are chrom, start and end,
# subsequent fields are optional: feature name, score, strand.
#
# We can get all of this info from the MIPGen output file. Chrom is the first part of the mip_key
# field, strand is the probe_strand (18th) field, mip_scan_start_position and mip_scan_stop_position
# (12th and 13th fields) give the start and stop positions and mip_name (field 20) gives the name.
# Score is aritrary
#
# *NBB* Bed files use 0-indexing

TARGETS_BED=$(addprefix $(VCF_DIR)/,targets.bed)

$(TARGETS_BED): $(CHOSEN_MIPS) | $(VCF_DIR)
	tail -n50 $(CHOSEN_MIPS) \
	| cut -f 1,12,13,18,20 \
	| cut -d "/" -f 1,3 \
	| sed s/\\/.// \
	| awk '{print $$1,$$2 - 1,$$3 - 1,$$5,100,$$4}' \
	| tr " " "\t" > $(TARGETS_BED)



#tail -n50 SoybeanLooper.Selected.MIPs.txt | cut -f 1,12,13,18,20 | cut -d "/" -f 1,3 | sed s/\\/.// | awk '{print $1,$1,$2,$3,$4,$5}' | cut -d ":" -f2,3 | sed 's/-/ /' | awk '{print $3,$4 - $1,($2 - $1) - ($2 - $5),$7,100,$6}' | tr " " "\t" > vcf/targets.bed
# NBB Tablet parses BED files as GFF3 and expects 9 columns

.PHONY: targets_bedfile

targets_bedfile: $(TARGETS_BED)

VCF_UNFILTERED=$(addprefix $(VCF_DIR)/,unfiltered.vcf)

# Call polymorphisms with freebayes
$(VCF_UNFILTERED): $(SORTED_ALIGNMENTS) $(MIPS_FASTA) $(TARGETS_BED) | $(VCF_DIR)
	freebayes -f $(MIPS_FASTA) -t $(TARGETS_BED) $(SORTED_ALIGNMENTS) > $(VCF_UNFILTERED)

.PHONY: call_polymorphisms

call_polymorphisms: $(VCF_UNFILTERED)

# GATK requires a sequence dictionary to do VCF filtering

GATK_DICT=$(subst .fasta,.dict,$(MIPS_FASTA))

$(GATK_DICT): $(MIPS_FASTA)
	gatk CreateSequenceDictionary -R $(MIPS_FASTA) -O $(GATK_DICT)

# Just the SNPS

VCF_SNPS=$(addprefix $(VCF_DIR)/,snps.vcf)

$(VCF_SNPS): $(VCF_UNFILTERED) $(GATK_DICT)
	gatk SelectVariants -R $(MIPS_FASTA) -V $(VCF_UNFILTERED) -O $(VCF_SNPS) -select-type SNP

# Filtered SNPs

MIN_POLYMORPHISM_QUAL=30

VCF_SNPS_FILTERED=$(addprefix $(VCF_DIR)/,filtered.snps.vcf)

$(VCF_SNPS_FILTERED): $(VCF_SNPS)
	gatk VariantFiltration \
	-R $(MIPS_FASTA) \
	-V $(VCF_SNPS) \
	-O $(VCF_SNPS_FILTERED) \
	--filter-expression "QUAL < 30" \
	--filter-name "qual" \
	--filter-expression "AB < 0.2 || AB > 0.8" \
	--filter-name "all-bal" \

# Filter for individuals with > 50% genotyped loci and loci where < 80% of individuals genotyped

ACCEPTED_INDS=$(addprefix $(VCF_DIR)/,accepted.inds)

$(ACCEPTED_INDS): $(VCF_SNPS_FILTERED)
	vcftools --vcf $(VCF_SNPS_FILTERED) --stdout --missing-indv --remove-filtered-all | \
	awk '$$5 < 0.5 {print $$1}' > $(ACCEPTED_INDS)

VCF_SNPS_INDS_FILTERED=$(addprefix $(VCF_DIR)/,filtered.snps.inds.vcf)

$(VCF_SNPS_INDS_FILTERED): $(VCF_SNPS_FILTERED) $(ACCEPTED_INDS)
	vcftools --vcf $(VCF_SNPS_FILTERED) \
	--stdout \
	--recode \
	--max-missing 0.8 \
	--remove-filtered-all \
	--keep $(ACCEPTED_INDS) \
	> $(VCF_SNPS_INDS_FILTERED)


VCF_HAPS=$(addprefix $(VCF_DIR)/,haplotypes.vcf)

# freebayes requires *input* VCFs to be bgzippes/tabixed

VCF_SNPS_INDS_FILTERED_GZ=$(addsuffix .gz, $(VCF_SNPS_INDS_FILTERED))

$(VCF_SNPS_INDS_FILTERED_GZ): $(VCF_SNPS_INDS_FILTERED)
	bgzip -f $(VCF_SNPS_INDS_FILTERED)
	tabix -p vcf $(VCF_SNPS_INDS_FILTERED_GZ)


$(VCF_HAPS): $(VCF_SNPS_INDS_FILTERED_GZ)
	freebayes -f $(MIPS_FASTA) --haplotype-length 100 --haplotype-basis-alleles $(VCF_SNPS_INDS_FILTERED_GZ) $(SORTED_ALIGNMENTS) > $(VCF_HAPS)

.PHONY: test

test: $(VCF_HAPS)
