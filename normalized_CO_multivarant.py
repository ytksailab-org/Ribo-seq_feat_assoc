import os
import sys
import csv
import numpy as np
import pandas as pd
import sys
import os.path
import scipy.stats as stats
import csv
import re
import os
from pathlib import Path
#python normalized_CO_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all out 2> log.txt


input_codon_wave = sys.argv[1]
input_path = sys.argv[2]

# maximum length of gene (aa)
MaxLen = 5000
# mean_co[i] will have the mean of CO at position i over all genes
mean_abs_co = [0] * MaxLen
mean_rel_co = [0] * MaxLen
# sum_co[i] will have the sum of CO at position i over all genes
sum_abs_co = [0] * MaxLen
sum_rel_co = [0] * MaxLen
# num_co[i] will have the number of genes having CO at position i
num_abs_co = [0] * MaxLen
num_rel_co = [0] * MaxLen
# Tips:
# x = [0] * L
# makes a list x whose length is L and whose values are all 0s



with open(input_codon_wave) as f:
    for line2 in f.readlines():
        liness = line2.split()
        geneid = liness[0]
        contactfile= Path("/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb.contact.order.fast.use".format(geneid))
        if not contactfile.is_file():
            continue
        contactfile = "/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb.contact.order.fast.use".format(geneid)
        abs_co = []
        rel_co = []

        with open(contactfile,"r") as contact:
            reader2 = csv.reader(contact, delimiter='\t')
            for row2 in reader2:
                abs_co.append(row2[1])
                rel_co.append(row2[2])
        if len(abs_co) > MaxLen:
            sys.exit(geneid + ' has CO longer than the maximum length of gene')
        for i in range(len(abs_co)):
            if not abs_co[i] == 'NA':
                sum_abs_co[i] += float(abs_co[i])
                num_abs_co[i] += 1
            if not rel_co[i] == 'NA':
                sum_rel_co[i] += float(rel_co[i])
                num_rel_co[i] += 1


# Make mean_co
for i in range(MaxLen):
    if num_abs_co[i] == 0:
        mean_abs_co[i] = 'NA'
    else:
        mean_abs_co[i] = sum_abs_co[i] / num_abs_co[i]
    if num_rel_co[i] == 0:
        mean_rel_co[i] = 'NA'
    else:
        mean_rel_co[i] = sum_rel_co[i] / num_rel_co[i]
for i in range(MaxLen):
    print('mean_co at ' + str(i) + '\t' + str(mean_abs_co[i]) + '\t' + str(mean_rel_co[i]), file=sys.stderr)


with open(input_codon_wave) as f:
    for line2 in f.readlines():
        liness = line2.split()
        geneid = liness[0]
        contactfile= Path("/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb.contact.order.fast.use".format(geneid))
        if not contactfile.is_file():
            continue
        contactfile = "/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb.contact.order.fast.use".format(geneid)
        abs_co = []
        rel_co = []
        with open(contactfile, "r") as contact1:
            reader2 = csv.reader(contact1, delimiter='\t')
            for row2 in reader2:
                abs_co.append(row2[1])
                rel_co.append(row2[2])
        if len(abs_co) > MaxLen:
            sys.exit(geneid + ' has CO longer than the maximum length of gene')
        norm_abs_co = [0] * len(abs_co)
        norm_rel_co = [0] * len(rel_co)
        for i in range(len(abs_co)):
            if abs_co[i] == 'NA':
                norm_abs_co[i] = 'NA'
            elif mean_abs_co[i] == 'NA':
                sys.exit('mean_co is not available for position ' + str(i))
            elif mean_abs_co[i] == 0.0 and float(abs_co[i]) == 0.0:
                norm_abs_co[i] = str(1.0)
            elif mean_abs_co[i] == 0.0 and float(abs_co[i]) != 0.0:
                sys.exit('illegal division by zero at ' + geneid + ' position ' + str(i))
            else:
                norm_abs_co[i] = str(float(abs_co[i]) / mean_abs_co[i])
            if rel_co[i] == 'NA':
                norm_rel_co[i] = 'NA'
            elif mean_rel_co[i] == 'NA':
                sys.exit('mean_co is not available for position ' + str(i))
            elif mean_rel_co[i] == 0.0 and float(rel_co[i]) == 0.0:
                norm_rel_co[i] = str(1.0)
            elif mean_rel_co[i] == 0.0 and float(rel_co[i]) != 0.0:
                sys.exit('illegal division by zero at ' + geneid + ' position ' + str(i))
            else:
                norm_rel_co[i] = str(float(rel_co[i]) / mean_rel_co[i])
        af_dir = '/home/bbian/Data_all_calibrate/result/CO_normalization/' + input_path
        output_co = af_dir + '/' + geneid + '.norm.txt'
        with open(output_co, 'w') as f3:
            writer = csv.writer(f3, delimiter='\t')
            for i in range(len(norm_abs_co)):
                writer.writerow([str(norm_abs_co[i]), str(norm_rel_co[i])])

