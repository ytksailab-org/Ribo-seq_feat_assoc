import sys
import csv
import re
import os
import scipy.stats as stats
from pathlib import Path
import os.path
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


#python collect_ASA_single_genes.py normalized.byRNA.Riboseq.codon.wave


input_data_wave = sys.argv[1]
input_path = sys.argv[2]
#output_extract = sys.argv[2]


with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        amino = line[1].split(',')
        reads = line[2].split(',')

        structure_1 = []
        structure_3 = []
        structure_4 = []

        pdbfile1 = Path("/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb".format(geneid))
        if not pdbfile1.is_file():
            continue
        pdbfile = "/Your/work/path/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb".format(geneid)
        p = PDBParser()
        structure = p.get_structure("pdbfile", pdbfile)
        model = structure[0]
        dssp = DSSP(model, pdbfile, dssp='~/local/conda/anaconda3-220108/bin/mkdssp')
        mylist = list(dssp.keys())
        for i in range(len(mylist)):
            a_key = mylist[i]
            outdata = dssp[a_key]
            structure_1.append(outdata[0])
            structure_3.append(outdata[2])
            structure_4.append(outdata[3])
        print(len(structure_1),flush=True)
        print(len(amino),flush=True)
        print(geneid,flush=True)
#    with open(output_extract, mode='w') as output:
        if len(structure_1)==len(amino)-1:
                for i in range(len(amino)-1):
                    daf = ("{}\t{}".format(structure_4[i],reads[i],))
                    print(daf,file=open("/Your/work/path/result/extract_ASA_each_gene_correlation/"+input_path+"/%s.txt" % geneid, "a"))
        else:
            continue
