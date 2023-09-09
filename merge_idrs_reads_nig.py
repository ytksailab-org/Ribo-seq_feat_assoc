import sys
import csv
from pathlib import Path



input_wave = sys.argv[1]
input_path = sys.argv[2]

with open(input_wave) as f:
    for line in f.readlines():
        lines = line.split()
        uniprotid= lines[0]
        count_reads= lines[2].split(',')
        #print(len(count_reads))
        idrs=[]
        idrfile = Path("/Your/work/path/result/IDRs_calculate/iupred2a/"+input_path+"/AF-{}-F1-model_v2.fa.txt".format(uniprotid))
        if not idrfile.is_file():
            continue
        idrfile = "/Your/work/path/result/IDRs_calculate/iupred2a/"+input_path+"/AF-{}-F1-model_v2.fa.txt".format(uniprotid)
        with open (idrfile, "r") as idr:
            for liness in idr.readlines():
                if liness[0] == "#":
                    continue
                line2 = liness.split()
                idrs.append(line2[2])
        print(len(idrs))
        print(len(count_reads))
        print(uniprotid)
        if len(idrs)== len(count_reads)-1:
            for i in range (len(count_reads)-1):
                daf = ("{}\t{}".format(idrs[i], count_reads[i], ))
                print(daf, file=open("/Your/work/path/result/merge_idrs_reads/"+input_path+"/%s.txt" % uniprotid,"a"))
            else:
                continue
