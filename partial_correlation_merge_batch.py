import os
import sys
import csv
import numpy as np
#import pandas as pd
import sys
import os.path
#import scipy.stats as stats
import csv
import re
import os
from pathlib import Path
# python normalized_CO_counts_merge.py normalized_genename_yeast_K15_Riboseq.codon.wave.txt

input_data_wave = sys.argv[1]
input_path = sys.argv[2]

#data2=[]
#input_data = sys.argv[1]
with open(input_data_wave) as f:
    for line in f.readlines():
        line1 = line.split()
        geneid = line1[0]
        count_reads = line1[2].split(',')
        #my_file = Path("/home/bbian/Data_all_calibrate/result/CO_normalization/"+input_path+"/{}.norm.txt".format(geneid))
        my_file1 = Path("/home/bbian/Data_all_calibrate/result/extract_ASA_each_gene_correlation/"+input_path+"/{}.txt".format(geneid))
        my_file2 = Path("/home/bbian/Data_all_calibrate/result/merge_idrs_reads/"+input_path+"/{}.txt".format(geneid))
        my_file3 = Path("/home/bbian/Data_all_calibrate/result/merge_reads_codonUsage/"+input_path+"/{}.txt".format(geneid))
        my_file4 = Path("/home/bbian/Data_all_calibrate/result/normalized_absolute_contact_order/"+input_path+"/{}.txt".format(geneid))
        my_file5 = Path("/home/bbian/Data_all_calibrate/result/relative_contact_order/"+input_path+"/{}.txt".format(geneid))
        my_file6 = Path("/home/bbian/Data_all_calibrate/result/absolute_contact_order/"+input_path+"/{}.txt".format(geneid))

        #my_file = Path("/Users/bianbian/Desktop/riboseq_tools/translation/normalization_counts_merge/{}.norm.txt".format(geneid))

        if not my_file1.is_file():
            continue
        ASA = []
        with open(my_file1, "r") as t:
            for line1s in t.readlines():
                line3 = line1s.split()
                ASA.append(line3[0])
        if not my_file1.is_file():
            continue

        idr = []
        with open(my_file2, "r") as t1:
            for line1s2 in t1.readlines():
                line4 = line1s2.split()
                idr.append(line4[0])
        if not my_file2.is_file():
            continue

        codonUsage = []
        with open(my_file3, "r") as t2:
            for line1s3 in t2.readlines():
                line5 = line1s3.split()
                codonUsage.append(line5[1])
        if not my_file3.is_file():
            continue

        normalized_CO = []
        with open(my_file4, "r") as t3:
            for line1s4 in t3.readlines():
                line6 = line1s4.split()
                normalized_CO.append(line6[0])
        if not my_file4.is_file():
            continue

        relative_CO = []
        with open(my_file5, "r") as t4:
            for line1s4 in t4.readlines():
                line7 = line1s4.split()
                relative_CO.append(line7[0])
        if not my_file5.is_file():
            continue

        absolute_CO = []
        with open(my_file6, "r") as t5:
            for line1s5 in t5.readlines():
                line8 = line1s5.split()
                absolute_CO.append(line8[0])
        if not my_file6.is_file():
            continue

        print(geneid)
        print(len(count_reads))
        print(len(ASA))
        print(len(idr))
        print(len(codonUsage))
        print(len(normalized_CO))
        print(len(relative_CO))
        print(len(absolute_CO))
        if len(ASA) == len(idr) == len(codonUsage)-1 == len(count_reads) - 1==len(normalized_CO)==len(relative_CO)==len(absolute_CO):
            for i in range(len(count_reads) - 1):
                daf = ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(count_reads[i], ASA[i], idr[i], codonUsage[i], normalized_CO[i], relative_CO[i], absolute_CO[i], ))
                print(daf, file=open("/home/bbian/Data_all_calibrate/result/partial_test/"+input_path+"/%s.txt" % geneid, "a"))
        else:
            continue
