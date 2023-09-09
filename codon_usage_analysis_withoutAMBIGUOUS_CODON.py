import sys
import os
import csv
import re


input_codon_usage_table = sys.argv[1]  # codon usage table
input_codon_wave = sys.argv[2]  #  codon waves
input_path = sys.argv[3]

codon2frequency = {}

with open(input_codon_usage_table) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        codon = row[0]
        frequency = row[1]
        codon2frequency[codon] = frequency
    #print(codon2frequency)

with open(input_codon_wave) as f:
    for line2 in f.readlines():
        liness = line2.split()
        uniprotid = liness[0]
        #print(uniprotid)
        count_reads = liness[2].split(',')
        amino = liness[1].split(',')
        #print(count_reads)
        codons= liness[3].split(',')
        print(codons)
        print(count_reads)
        codon_upper = []
        for x in codons:
            codon_upper.append(x.upper())
        print(codon_upper)
        codon_frequency = []
        for i in codon_upper:
            if i in codon2frequency:
                codon_frequency.append(codon2frequency[i])
            else:
                codon_frequency.append('AMBIGUOUS_CODON')
        print(codon_frequency)
        print(len(amino))
        print(len(count_reads))
        print(len(codon_frequency))

        if len(count_reads) == len(codon_frequency):
            for i in range(len(count_reads)):
                if codon_frequency[i] != 'AMBIGUOUS_CODON':
                    reads_frequency = ("{}\t{}".format(count_reads[i], codon_frequency[i], ))
                print(reads_frequency, file=open("/home/bbian/Data_all_calibrate/result/merge_reads_codonUsage/"+input_path+"/%s.txt" % uniprotid, "a"))
        else:
            continue







