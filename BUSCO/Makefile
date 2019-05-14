##################
# BUSCO analysis #
##################

SoybeanLooperScaffolds.fasta: ../genome/SoybeanLooperScaffolds.fasta
	ln -f ../genome/SoybeanLooperScaffolds.fasta ./

BUSCO_OUT_BASE=busco_looper

LINEAGE_DIR=endopterygota_odb9

$(LINEAGE_DIR):
	wget https://busco.ezlab.org/datasets/endopterygota_odb9.tar.gz;
	tar xzvf endopterygota_odb9.tar.gz

export AUGUSTUS_CONFIG_PATH=/home/nick/anaconda3/envs/SoybeanLooper/config

# BUSCO produces a number of output files, pick one as the target

BUSCO_OUT=$(addprefix $(BUSCO_DIR)/,$(subst XXXX,$(BUSCO_OUT_BASE),short_summary_XXXX.txt))

$(BUSCO_OUT): SoybeanLooperScaffolds.fasta $(LINEAGE_DIR)
	run_BUSCO.py -i SoybeanLooperScaffolds.fasta -o $(BUSCO_OUT_BASE) -l $(LINEAGE_DIR) -m genome -c 24


.PHONY: busco

busco: $(BUSCO_OUT)
