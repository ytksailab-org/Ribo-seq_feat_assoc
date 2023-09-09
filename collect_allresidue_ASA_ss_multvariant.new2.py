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


#python collect_allresidue_ASA.py normalization.codon.wave collect_all_residue_ASA_use.txt


input_data_wave = sys.argv[1]
input_path = sys.argv[2]
output_extract = sys.argv[3]


with open(input_data_wave) as f:
    for line in f.readlines():
        line = line.split()
        geneid = line[0]
        amino = line[1].split(',')
        reads = line[2].split(',')
        print(geneid)
        print(len(reads))
        print(len(amino))

        structure_1 = []
        structure_3 = []
        structure_4 = []

        pdbfile1 = Path("/home/bbian/Data_all/raw_data/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb".format(geneid))
        if not pdbfile1.is_file():
            continue
        pdbfile = "/home/bbian/Data_all/raw_data/Alpha-fold/"+input_path+"/AF-{}-F1-model_v2.pdb".format(geneid)
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
        print(len(structure_1))
        print(len(amino))
        print(len(reads))
        print(geneid)


        if len(structure_1)==len(amino)-1:
            for i in range(len(amino)-1):
                print("{}\t{}\t{}\t{}\t{}".format(structure_1[i],amino[i],structure_3[i],structure_4[i],reads[i],),file= open(output_extract,"a"))
        else:
            continue

