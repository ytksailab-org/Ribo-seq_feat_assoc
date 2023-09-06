import sys
import csv
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

#python filter_low_reads_multivarant.py multi.variant.codon.wave.adjusted.format filtered.multi.variant.codon.wave.adjusted.format

input_data_wave = sys.argv[1]

output_new_wave = sys.argv[2]

with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        amino_seq= line[1].split(',')
        reads = line[2].split(',')

        if len(amino_seq)==len(reads):
            print("{}\t{}\t{}\t{}".format(line[0], line[1], line[2], line[3], ), file=open(output_new_wave, "a"))
        else:
            print(geneid)
            continue
