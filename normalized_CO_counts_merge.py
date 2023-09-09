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


input_data_wave = sys.argv[1]
input_path = sys.argv[2]


with open(input_data_wave) as f:
    for line in f.readlines():
        line1 = line.split()
        geneid = line1[0]
        count_reads = line1[2].split(',')
        my_file = Path("/Your/work/path/result/CO_normalization/"+input_path+"/{}.norm.txt".format(geneid))

        if not my_file.is_file():
            continue
        contact_order = []
        with open(my_file, "r") as t:
            for line1s in t.readlines():
                line3 = line1s.split()
                contact_order.append(line3[0])

        print(len(count_reads))
        print(len(contact_order))
        print(geneid)
        if len(contact_order) == len(count_reads) - 1:
            for i in range(len(count_reads) - 1):
                daf = ("{}\t{}".format(contact_order[i], count_reads[i], ))
                print(daf, file=open("/Your/work/path/result/normalized_absolute_contact_order/"+input_path+"/%s.txt" % geneid, "a"))
        else:
            continue
