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

GENOME_RAW_FASTQ=$(addprefix GenomeRawData/, MA22_S1_L001_R1_001.fastq.gz  MA22_S1_L001_R2_001.fastq.gz)

#
# Initial QC of the fastq files
#

QC_OUT_DIR=fastqc_out

FASTQC_OPTS=-o $(QC_OUT_DIR) -t 16

$(QC_OUT_DIR):
	if [ ! -d $(QC_OUT_DIR) ]; then mkdir $(QC_OUT_DIR); fi

fastqc_raw: $(QC_OUT_DIR)
	fastqc $(FASTQC_OPTS) $(GENOME_RAW_FASTQ)

test:
	echo $(GENOME_RAW_FASTQ)
