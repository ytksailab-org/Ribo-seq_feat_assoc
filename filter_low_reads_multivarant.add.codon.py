import sys
import csv
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

#python filter_low_reads_multivarant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60

input_data_wave = sys.argv[1]
#output_zero_name = sys.argv[2]
output_new_wave = sys.argv[2]
percent_high_zero = 0
percent_low_zero = 0

#with open(output_new_wave, mode='w') as fout:
    #writer = csv.writer(fout, delimiter='\t', lineterminator='\n')

with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        reads = line[2].split(',')
        counts = 0
        highzero_namelist = []
        genome_id = []
        gene_id =[]
        position = []
        readsnew = []
        codon = []
        amino_seq = []

        print(geneid)
        print (reads)
        print(len(reads))
        for i in reads:
            if (float(i) == 0):
                counts += 1
        print(counts)
        percent = counts/len(reads)
        print(percent)
        if percent >= 0.4:
            percent_high_zero +=1
            highzero_namelist.append(geneid)
            #print(highzero_namelist, file=open(output_zero_name, "a"))
            #outdata = [highzero_namelist]
            #writer.writerow(highzero_namelist)
            #with open(output_zero_name, mode='w') as output:
                #output.write(str(highzero_namelist))
            #print(highzero_namelist, file=open(output_zero_name, "a"))
        else:
            percent_low_zero += 1
            #genome_id.append(line[0])
            gene_id.append(line[0])
            #position.append(line[2])
            readsnew.append(line[2])
            #codon.append(line[4])
            amino_seq.append(line[1])
            outdata = [line[0],
                       line[1],
                       line[2],]
            #print(outdata, file=open(output_new_wave, "a"))
            print("{}\t{}\t{}\t{}".format(line[0], line[1], line[2], line[3],),file=open(output_new_wave, "a"))




    print(percent_high_zero,percent_low_zero)





