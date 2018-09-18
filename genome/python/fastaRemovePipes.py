# Ad-hoc script to remove pipe characters from identifiers in a fasta file
#
# Requires Biopython
#
# As a bonus, converts fasta files to 2-line format

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

infileName = sys.argv[1]

inFile = SeqIO.parse(infileName, "fasta")

for record in inFile:
    SeqIO.write(SeqRecord(record.seq, id = record.id.split("|")[0], description = ""), sys.stdout, "fasta-2line")


