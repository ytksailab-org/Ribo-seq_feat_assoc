import numpy as np
import pandas as pd
import sys
import os.path
import scipy.stats as stats
import csv
import re
import os
from pathlib import Path



input_data_wave = sys.argv[1]
input_path = sys.argv[2]


with open(input_data_wave) as f:
    for line2 in f.readlines():
        liness = line2.split()
        geneid = liness[0]
        count_reads = liness[2].split(',')

        #contact_order = []

        contactfile= Path("/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb.contact.order.fast.use".format(geneid))
        if not contactfile.is_file():
            continue
        contactfile = "/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb.contact.order.fast.use".format(geneid)
        contact_order = []
        
        with open(contactfile,"r") as contact:
            for line1s in contact.readlines():
                line3 = line1s.split()
                #output = line[0] + "\t" + line[1] + "\t" + line[2] + "\n"
                contact_order.append(line3[2])

        print(len(count_reads))
        print(len(contact_order))
        print(geneid)
        if len(contact_order) == len(count_reads) - 1:
            for i in range(len(count_reads) - 1):
                daf = ("{}\t{}".format(contact_order[i], count_reads[i], ))
                print(daf, file=open("/Your/work/path/result/relative_contact_order/"+input_path+"/%s.txt" % geneid, "a"))
        else:
            continue
