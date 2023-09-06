import sys
import csv
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

#python python scale_Riboseq_multivariant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60 yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all

input_data_wave = sys.argv[1]
#output_zero_name = sys.argv[2]
output_new_wave = sys.argv[2]


#with open(output_new_wave, mode='w') as fout:
    #writer = csv.writer(fout, delimiter='\t', lineterminator='\n')

with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        reads = line[2].split(',')
        counts_each_gene = 0
        readsnew = []


        print(geneid)
        print (reads)
        print(line[2])
        print(len(reads))
        for i in reads:
            counts_each_gene +=int(i)
        print(counts_each_gene)
        mean_counts = counts_each_gene/len(reads)
        print(mean_counts)
        for a in reads:
            new_counts = int(a) / mean_counts
            readsnew.append(new_counts)
        print(*readsnew)
        print(readsnew)

        read_use=""
        for i in readsnew:
            read_use =read_use + ","+str(i)
            Read_use = read_use[1:len(read_use)]
        print(Read_use)
        #print((line[0], line[1], line[2], readsnew, line[4], line[5],),file=open(output_new_wave, "a")
        print("{}\t{}\t{}\t{}".format(line[0], line[1], Read_use, line[3],),file=open(output_new_wave, "a"))
