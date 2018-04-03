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
