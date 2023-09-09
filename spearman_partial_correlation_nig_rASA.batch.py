#import numpy as np
import pandas as pd
import pingouin as pg
import sys
import csv
import re
import os

input_data_wave = sys.argv[1]
input_path = sys.argv[2]
#data2=[]
#input_data = sys.argv[1]
with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        try:
            
            with open("/Your/work/path/result/partial_test/"+input_path+"/{}.txt".format(geneid), "r") as f:
                a = pd.read_csv(f,sep='\t',header =None, names= ["ReadsCount","RelativeASA", "idr", "codon_usage", "normalized_CO", "relative_CO", "absolute_CO"])
                a1 = a.dropna()
                data1= pg.partial_corr(data=a1, x='ReadsCount', y='RelativeASA', covar=['idr', 'codon_usage', 'normalized_CO', 'relative_CO','absolute_CO'],method = 'spearman').round(3)
                shel= data1.iat[0, 1]
                #data2.append(data1.iat[0,1])
                #print(data1)
                print(shel)
                #print(data1.iat[0, 1])
                #print(geneid,data2)
            daf = ("{}\t{}".format(geneid, shel, ))
            print(daf, file=open("/Your/work/path/result/partial_test/spearman_partial/"+input_path+".partial.single_gene_corre_ASA.txt", "a"))
        except IOError:
            continue


