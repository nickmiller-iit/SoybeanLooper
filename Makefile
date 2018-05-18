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
# dir
#

REDUNDANSFILLEDSCAFFS=$(addsuffix /scaffolds.filled.fa, $(REDUNDANS_OUT))

SoybeanLooperScaffolds.fasta: $(REDUNDANSFILLEDSCAFFS)
	cp -L $(REDUNDANSFILLEDSCAFFS) SoybeanLooperScaffolds.fasta

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
	ln SoybeanLooperScaffolds.fasta $(BWAREF)


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
