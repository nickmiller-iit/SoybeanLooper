########################################################################################
#                                                                                      #
# Soybean looper analysis.                                                             #
#                                                                                      #
# Note that these analyses should be run under the "SoybeanLooper" conda environment   #
#                                                                                      #
########################################################################################


##############################
# Assembling genomic contigs #
##############################

GENOME_RAW_DIR=GenomeRawData
GENOME_RAW1=MA22_S1_L001_R1_001.fastq.gz
GENOME_RAW2=MA22_S1_L001_R2_001.fastq.gz
GENOME_RAW_ALL=$(addprefix $(GENOME_RAW_DIR)/, $(GENOME_RAW1) $(GENOME_RAW2))

#
# Initial QC of the fastq files
#

QC_OUT_DIR=fastqc_out

FASTQC_OPTS=-o $(QC_OUT_DIR) -t 16

$(QC_OUT_DIR):
	if [ ! -d $(QC_OUT_DIR) ]; then mkdir $(QC_OUT_DIR); fi

fastqc_raw: $(QC_OUT_DIR)
	fastqc $(FASTQC_OPTS) $(GENOME_RAW_ALL)

#
# Quality clipping and adaptor trimming of raw data
#

GENOME_TRIMMED_DIR=GenomeTrimmedData

$(GENOME_TRIMMED_DIR):
	if [ ! -d $(GENOME_TRIMMED_DIR) ]; then mkdir $(GENOME_TRIMMED_DIR); fi

GENOME_TRIMMED_1_BASE=$(addprefix $(GENOME_TRIMMED_DIR)/, $(GENOME_RAW1))

GENOME_TRIMMED_2_BASE=$(addprefix $(GENOME_TRIMMED_DIR)/, $(GENOME_RAW2))

GENOME_TRIMMED_1P=$(subst R1_001,R1_P,$(GENOME_TRIMMED_1_BASE))

GENOME_TRIMMED_1UP=$(subst R1_001,R1_UP,$(GENOME_TRIMMED_1_BASE)) 

GENOME_TRIMMED_2P=$(subst R2_001,R2_P,$(GENOME_TRIMMED_2_BASE))

GENOME_TRIMMED_2UP=$(subst R2_001,R2_UP,$(GENOME_TRIMMED_2_BASE)) 

GENOME_TRIMMED_ALL=$(GENOME_TRIMMED_1P) $(GENOME_TRIMMED_1UP) $(GENOME_TRIMMED_2P) $(GENOME_TRIMMED_2UP)

GENOME_TRIMS=ILLUMINACLIP:/home/nick/anaconda3/envs/SoybeanLooper/share/trimmomatic-0.36-5/adapters/NexteraPE-PE.fa:3:25:7 SLIDINGWINDOW:4:20 MINLEN:40

$(GENOME_TRIMMED_ALL): $(GENOME_TRIMMED_DIR) $(GENOME_RAW_ALL)
	trimmomatic PE -threads 16 $(GENOME_RAW_ALL) $(GENOME_TRIMMED_ALL) $(GENOME_TRIMS)

trim_genome: $(GENOME_TRIMMED_ALL)

#
# QC of trimmed data
#

fastqc_trimmed: $(QC_OUT_DIR) $(GENOME_TRIMMED_ALL)
	fastqc $(FASTQC_OPTS) $(GENOME_TRIMMED_ALL)

#
# Assemble with SPAdes
#

# SPAdes uses an output dir, so just set that as the target.

spades.assembly: $(GENOME_TRIMMED_ALL)
	spades.py -o spades.assembly -1 $(GENOME_TRIMMED_1P) -2 $(GENOME_TRIMMED_2P)

#
# Reduce redundancy and scaffold with redundans
#

# Redundans is not available on bioconda but it is available as a docker image
# Using the docker image to save installing a bunch of dependencies including python 2

#
# We will need to give docker access to the contigs from SPAdes and the trimmed fastq files
# after some playing around, it looks like docker is a bit ticky about mounting multiple different
# volumes. Solution is to set up a redundans dir, and link in the required files
#

REDUNDANS_DIR=redundans

$(REDUNDANS_DIR):
	if [ ! -d $(REDUNDANS_DIR) ]; then mkdir $(REDUNDANS_DIR); fi

REDUNDANS_FASTA_IN=$(REDUNDANS_DIR)/contigs.fasta

$(REDUNDANS_FASTA_IN): $(REDUNDANS_DIR) spades.assembly/contigs.fasta
	ln spades.assembly/contigs.fasta $(REDUNDANS_FASTA_IN)

REDUNDANS_TRIMMED_1P=$(subst $(GENOME_TRIMMED_DIR),$(REDUNDANS_DIR),$(GENOME_TRIMMED_1P))

$(REDUNDANS_TRIMMED_1P): $(REDUNDANS_DIR) $(GENOME_TRIMMED_1P)
	ln $(GENOME_TRIMMED_1P) $(REDUNDANS_TRIMMED_1P)

REDUNDANS_TRIMMED_2P=$(subst $(GENOME_TRIMMED_DIR),$(REDUNDANS_DIR),$(GENOME_TRIMMED_2P))

$(REDUNDANS_TRIMMED_2P): $(REDUNDANS_DIR) $(GENOME_TRIMMED_2P)
	ln $(GENOME_TRIMMED_2P) $(REDUNDANS_TRIMMED_2P)

REDUNDANS_INPUTS=$(REDUNDANS_FASTA_IN) $(REDUNDANS_TRIMMED_1P) $(REDUNDANS_TRIMMED_2P)

DCKR_REDUNDANS_CMD=docker run -it -v $(CURDIR)/$(REDUNDANS_DIR):/working  lpryszcz/redundans /root/src/redundans/redundans.py

REDUNDANS_ARGS=-v -t 16 -f /working/contigs.fasta -i $(subst $(REDUNDANS_DIR),/working,$(REDUNDANS_TRIMMED_1P)) $(subst $(REDUNDANS_DIR),/working,$(REDUNDANS_TRIMMED_2P)) -o /working/out

REDUNDANS_OUT=redundans/out

$(REDUNDANS_OUT): $(REDUNDANS_INPUTS)
	$(DCKR_REDUNDANS_CMD) $(REDUNDANS_ARGS)

run_redundans: $(REDUNDANS_OUT)

#
# Redundans symlinks its output scaffolds. Make a copy of the filled scaffolds file in the top level working
# dir. While we are doing this, get read of pipe chars in fasta ID lines, which cause downstream headaches
# and convert records to 2-line fasta format.
#

REDUNDANSFILLEDSCAFFS=$(addsuffix /scaffolds.filled.fa, $(REDUNDANS_OUT))

SoybeanLooperScaffolds.fasta: $(REDUNDANSFILLEDSCAFFS)
	python python/fastaRemovePipes.py $(REDUNDANSFILLEDSCAFFS) > SoybeanLooperScaffolds.fasta
#	cp -L $(REDUNDANSFILLEDSCAFFS) SoybeanLooperScaffolds.fasta

#
# Mask repeats
#

#
# A 2-step process: First find repeats with RepeatModeller, then mask them with RepeatMasker
#

# Input is filled scaffolds from redundans

REDUNDANS_SCAFFS=$(addsuffix /scaffolds.filled.fa, $(REDUNDANS_OUT))

RPTS_DIR=repeatmask

# Copy and rename the redundands scaffolds file to something more meaningful

UNMASKED=$(addsuffix /SoybeanLooper.fasta, $(RPTS_DIR))

$(UNMASKED): $(REDUNDANS_SCAFFS)
	if [ ! -d $(RPTS_DIR) ]; then mkdir $(RPTS_DIR); fi
	cp $(REDUNDANS_SCAFFS) $(UNMASKED)


#
# RepeatModeler requires a blast database of the genome.
# the bundled BuildDatabase utility will do this.
# The blast db consists of multiple files in the form dbname.*
# We will pick one of these file as the target for make.
#

RPT_DB=$(addsuffix /SoybeanLooper, $(RPTS_DIR))

RPT_DB_TARGET=$(addsuffix .nin, $(RPT_DB))

$(RPT_DB_TARGET): $(UNMASKED)
	BuildDatabase -name $(RPT_DB) -engine ncbi $(UNMASKED)

#
# RepeatModeler also generates > 1 output file
# Use the consensus sequences ouput as the target
#

RPT_SEQS=$(addsuffix -families.fa, $(RPT_DB))

$(RPT_SEQS): $(RPT_DB_TARGET)
	RepeatModeler -engine ncbi -pa 16 -database $(RPT_DB)

repeat_models: $(RPT_SEQS)


#
# Align reads with bwa mem
#

BWADIR=bwa

#
# Set up db. Link in the assembled scaffolds
#

BWAREF=$(addsuffix /SoybeanLooperScaffolds.fasta, $(BWADIR))

$(BWAREF): SoybeanLooperScaffolds.fasta
	if [ ! -d $(BWADIR) ]; then mkdir $(BWADIR); fi
	ln -f SoybeanLooperScaffolds.fasta $(BWAREF)


# BWA index consists of several files with different extensions, use .bwt as a "trigger"


BWAIDX=$(addsuffix .bwt, $(BWAREF))

$(BWAIDX): $(BWAREF)
	bwa index $(BWAREF)

#
# Align with bwa mem
#

BWAALIGN=$(addsuffix .bam, $(BWAREF))

$(BWAALIGN): $(BWAIDX) $(GENOME_TRIMMED_1P) $(GENOME_TRIMMED_2P)
	bwa mem -t 24 $(BWAREF) $(GENOME_TRIMMED_1P) $(GENOME_TRIMMED_2P) | samtools view -b -o $(BWAALIGN)

bwa_alignment: $(BWAALIGN)

#
# Get rid of PCR duplicates. This is a bit of a business due to different sorting requirements
#

DEDUPALIGN=$(addsuffix .dedup.bam, $(BWAREF))

$(DEDUPALIGN): $(BWAALIGN)
	samtools sort -n -@ 24 $(BWAALIGN) > $(BWADIR)/tmp1.bam
	samtools fixmate -m $(BWADIR)/tmp1.bam $(BWADIR)/tmp2.bam
	samtools sort -@ 24 $(BWADIR)/tmp2.bam > $(BWADIR)/tmp3.bam
	samtools markdup -r $(BWADIR)/tmp3.bam $(DEDUPALIGN)
	rm $(BWADIR)/tmp*.bam

bwa_alignment_dedup: $(DEDUPALIGN)

#
# Extract read depth stats, including sites with 0 depth
#

DEDUP_DEPTH_STATS=$(addsuffix /depth.stats, $(BWADIR))

$(DEDUP_DEPTH_STATS): $(DEDUPALIGN)
	samtools depth -a $(DEDUPALIGN) > $(DEDUP_DEPTH_STATS)

depth_stats: $(DEDUP_DEPTH_STATS)

#
# Add an index so we can view the alignment with Tablet
#

DEDUPALIGN_IDX=$(addsuffix .dedup.bam.bai, $(BWAREF))

$(DEDUPALIGN_IDX): $(DEDUPALIGN)
	samtools index $(DEDUPALIGN) $(DEDUPALIGN_IDX)

BWA_REF_IDX=$(addsuffix .fai, $(BWAREF))

$(BWA_REF_IDX): $(BWAREF)
	samtools faidx $(BWAREF)

index_bwa_alignment_dedup: $(DEDUPALIGN_IDX) $(BWA_REF_IDX)


#######################################
# Identifying target regions for MIPs #
#######################################

BEDDIR=bed

$(BEDDIR):
	if [ ! -d $(BEDDIR) ]; then mkdir $(BEDDIR); fi

#
# Make a BED file of contigs >= 1 kb
#
# Because MIPgen will fail if we try and tile probes at the very start or end of the contif, we do not include
# the first  and last 200 bp of the contig in the BED file. For example, a 1000 bp contig will yeild a line in
# th BED file indicating the middle 600bp
#
1KBCTGS_BED=$(addsuffix /1kbContigs.bed, $(BEDDIR))

$(1KBCTGS_BED): $(DEDUPALIGN_IDX) | $(BEDDIR)
	samtools idxstats $(DEDUPALIGN) | awk '$$2 >= 1000' | awk '{OFS = "\t"; print $$1, 199, $$2-201}' > $(1KBCTGS_BED)

SINGLECOPY_BED=$(addsuffix /singlecopy.bed, $(BEDDIR))

#
# Make a bedfile of regions > 200bp where coverage is <= 50 reads pre site
#

$(SINGLECOPY_BED): $(DEDUPALIGN_IDX) | $(BEDDIR)
	bedtools genomecov -bg -ibam $(DEDUPALIGN) | awk '$$4 <= 50' | cut -f 1-3 | bedtools merge | awk '$$3 - $$2 > 199' > $(SINGLECOPY_BED)

#
# Make a bedfile of the intersection of "single copy" regions and 1kb contigs
#

SINGLECOPY_1KBCTGS_BED=$(addsuffix /1kbctgs.singlecopy.bed, $(BEDDIR))

$(SINGLECOPY_1KBCTGS_BED): $(1KBCTGS_BED) $(SINGLECOPY_BED)
	bedtools intersect -a $(1KBCTGS_BED) -b $(SINGLECOPY_BED) > $(SINGLECOPY_1KBCTGS_BED)

singlecopy_bedfile: $(SINGLECOPY_1KBCTGS_BED)


#########################################
# Generating MIPs sequences with mipgen #
#########################################

#
# Download and build mipgen.
# Mipgen is not available via bioconda or (afaik) other pre-built binary sources
# Instead, download an build locally.
#
# We have to patch miopgen.cpp to avoid the normal behaviour of looking up sequence +/- 1 kb of the target
#

MIPGEN_BINDIR=MIPGEN

MIPGEN_BIN=$(addsuffix /mipgen, $(MIPGEN_BINDIR))

$(MIPGEN_BINDIR):
	git clone https://github.com/shendurelab/MIPGEN.git


$(MIPGEN_BIN): | $(MIPGEN_BINDIR)
	cd $(MIPGEN_BINDIR) &&\
	cp ../patches/mipgen.patch ./ &&\
	git apply mipgen.patch &&\
	$(MAKE)
#	$(MAKE) -C $(MIPGEN_BINDIR)

mipgen: $(MIPGEN_BIN)


#
# Make a BED file of targets for MIPGEN
#
# The actual set of targets is a sample of all possible targets. If we use everythin, the output
# files from MIPgen get up to several terabytes. The sample is 10000 stargets, chosen at random
# Using a small R script. Give the script a fixed seed to ensure repreatable picking of the same subset.
#

MIPGEN_DIR=mipgen

$(MIPGEN_DIR):
	if [ ! -d $(MIPGEN_DIR) ]; then mkdir $(MIPGEN_DIR); fi

MIPGEN_TARGETS=$(addsuffix /targets.bed, $(MIPGEN_DIR))


$(MIPGEN_TARGETS): $(SINGLECOPY_1KBCTGS_BED) | $(MIPGEN_DIR)
	bedtools makewindows -b $(SINGLECOPY_1KBCTGS_BED) -w 1 -s 500  > $(MIPGEN_DIR)/targets_all.bed
	Rscript ./R/subsampleBED.R $(MIPGEN_DIR)/targets_all.bed 10000 42 > $(MIPGEN_TARGETS)


#
# Make a Samtools/BWA-Indexed genome
#

MIPGEN_GENOME=$(addsuffix /SoybeanLooperScaffolds.fasta, $(MIPGEN_DIR))

$(MIPGEN_GENOME): SoybeanLooperScaffolds.fasta | $(MIPGEN_DIR)
	ln -f SoybeanLooperScaffolds.fasta $(MIPGEN_GENOME)

MIPGEN_GENOME_BWAIDX=$(addsuffix .bwt, $(MIPGEN_GENOME))

$(MIPGEN_GENOME_BWAIDX): $(MIPGEN_GENOME)
	bwa index $(MIPGEN_GENOME)

MIPGEN_GENOME_FAIDX=$(addsuffix .fai, $(MIPGEN_GENOME))

$(MIPGEN_GENOME_FAIDX): $(MIPGEN_GENOME)
	samtools faidx $(MIPGEN_GENOME)

#
# run MIPgen
#

MIPGEN_PROJ=$(addsuffix /SoybeanLooper, $(MIPGEN_DIR))

MIPGEN_OPTS=-min_capture_size 154 -max_capture_size 184 -tag_sizes 0,0 -bwa_threads 16

MIPGEN_OUT=$(addsuffix /$(addsuffix .picked_mips.txt, $(MIPGEN_PROJ)), $(MIPGEN_DIR))

$(MIPGEN_OUT): $(MIPGEN_BIN) $(MIPGEN_TARGETS) $(MIPGEN_GENOME_BWAIDX) $(MIPGEN_GENOME_FAIDX)
	$(MIPGEN_BIN) -regions_to_scan $(MIPGEN_TARGETS) -project_name $(MIPGEN_PROJ) $(MIPGEN_OPTS) -bwa_genome_index $(MIPGEN_GENOME)

mips_all: $(MIPGEN_OUT)

test: $(MIPGEN_TARGETS)
