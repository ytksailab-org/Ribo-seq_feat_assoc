import sys
import os
import csv
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

input_codon_wave = sys.argv[1]  # Human_K15_SRR1562539_codon.wave
input_geneid2genename = sys.argv[2]  # geneID2genename.table.txt
input_genename2uniprotid = sys.argv[3]  # human_genename2uniprotid.table.txt
AlphafoldPH = sys.argv[4]


# Algorithm:
# 1. Scan input_codon_wave to make following dictionaries.
#    A. From aaseq to a list of geneids (to know geneids with the same aaseq)
#    B. From geneid to codon wave positions
#    C. From geneid to codon wave counts
#    D. From geneid to codon wave codons
#    E. From geneid to aaseq
# 2. Make a dictionary from aaseq to codon wave counts.
#    For each aaseq in dictionary 1A,
#    add (or average) codon wave counts of corresponding geneids obtained from dictionary 1C.
# 3. Scan input_geneid2genename and input_genename2uniprotid to make following dictionaries:
#    A. From uniprotid to a list of genenames
#    B. From genename to a list of geneids
# 4. Using dictionaries 3A and 3B,
#    make a dictionary from uniprotid to a list of geneids
# 5. For each uniprotid in dictionary 4,
#    5.1 Obtain multiple aaseqs of corresponding geneids from dictionary 1E
#    5.2 Load aaseq from the uniprodid's AlphaFold pdb file
#    5.3 Compare each of aaseqs from 5.1 and aaseq from 5.2,
#    5.4 Select the most similar aaseq from 5.1
#    5.5 make a dictionary from uniprodid to the selected aaseq
# Then, output the results with your favorite format.
# For example, for each uniprotid, output followings things as one line
# - uniprotid
# - aaseq from dictionary 5.5 (not from pdb file)
# - multiple geneids (comma separated) from dictionary 1A (not from dictionary 4)
# - codon wave counts from dictionary 2
# - etc

csv.field_size_limit(1000000)

# 1. Scan input_codon_wave to make following dictionaries.
print('Starting Step 1', file=sys.stderr, flush=True)
#    A. From aaseq to a list of geneids (to know geneids with the same aaseq)
aaseq2geneids = {}
#    B. From geneid to codon wave positions
geneid2codon_wave_positions = {}
#    C. From geneid to codon wave counts
geneid2codon_wave_counts = {}
#    D. From geneid to codon wave codons
geneid2codon_wave_codons = {}
#    E. From from geneid to aaseq
geneid2aaseq = {}
with open(input_codon_wave) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        geneid = row[1]
        codon_wave_positions = row[2].split(',')
        codon_wave_positions = [int(s) for s in codon_wave_positions]
        codon_wave_counts = row[3].split(',')
        codon_wave_counts = [int(s) for s in codon_wave_counts]
        codon_wave_codons = row[4].split(',')
        codon_wave_aas = row[5].split(',')
        aaseq = "".join(codon_wave_aas)
        if aaseq not in aaseq2geneids:
            aaseq2geneids[aaseq] = []
        aaseq2geneids[aaseq].append(geneid)
        geneid2codon_wave_positions[geneid] = codon_wave_positions
        geneid2codon_wave_counts[geneid] = codon_wave_counts
        geneid2codon_wave_codons[geneid] = codon_wave_codons
        geneid2aaseq[geneid] = aaseq

# 2. Make a dictionary from aaseq to codon wave counts.
print('Starting Step 2', file=sys.stderr, flush=True)
aaseq2codon_wave_counts = {}
#    For each aaseq in dictionary 1A,
#    add (or average) codon wave counts of corresponding geneids obtained from dictionary 1C.
for aaseq in aaseq2geneids.keys():
    wave_len = -1
    codon_wave_counts = []
    for geneid in aaseq2geneids[aaseq]:
        if wave_len == -1:
            wave_len = len(geneid2codon_wave_counts[geneid])
            codon_wave_counts = [0] * wave_len
        if wave_len != len(geneid2codon_wave_counts[geneid]):
            sys.exit('inconsistent length of wave for ' + geneid)
        for i in range(wave_len):
            codon_wave_counts[i] += geneid2codon_wave_counts[geneid][i]
    aaseq2codon_wave_counts[aaseq] = codon_wave_counts

# 3. Scan input_geneid2genename and input_genename2uniprotid to make following dictionaries:
print('Starting Step 3', file=sys.stderr, flush=True)
#    A. From uniprotid to a list of genenames
uniprotid2genenames = {}
#    B. From genename to a list of geneids
genename2geneids = {}
with open(input_genename2uniprotid) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        genenames = row[0].split(',')
        uniprotids = row[1].split(',')
        for genename in genenames:
            for uniprotid in uniprotids:
                if uniprotid not in uniprotid2genenames:
                    uniprotid2genenames[uniprotid] = []
                uniprotid2genenames[uniprotid].append(genename)
with open(input_geneid2genename) as f:
    reader = csv.reader(f, delimiter='\t')
    #header = next(reader)
    for row in reader:
        geneids = row[0].split(',')
        genenames = row[1].split(',')
        for geneid in geneids:
            for genename in genenames:
                if genename not in genename2geneids:
                    genename2geneids[genename] = []
                genename2geneids[genename].append(geneid)

# 4. Using dictionaries 3A and 3B,
#    make a dictionary from uniprotid to a list of geneids
print('Starting Step 4', file=sys.stderr, flush=True)
uniprotid2geneids = {}
for uniprotid in uniprotid2genenames.keys():
    for genename in uniprotid2genenames[uniprotid]:
        for geneid in genename2geneids[genename]:
            if uniprotid not in uniprotid2geneids:
                # We use dictionary rather than list to prevent duplicate entries
                uniprotid2geneids[uniprotid] = {}
            uniprotid2geneids[uniprotid][geneid] = 1

print('Correspondene between uniprotid and geneid', file=sys.stderr, flush=True)
for uniprotid in uniprotid2geneids.keys():
    print('\t'.join([uniprotid, ','.join(uniprotid2geneids[uniprotid].keys())]), file=sys.stderr, flush=True)

# 5. For each uniprotid in dictionary 4,
print('Starting Step 5', file=sys.stderr, flush=True)
uniprotid2aaseq = {}
uniprotid2best_geneid = {}
for uniprotid in uniprotid2geneids:
    #    5.1 Obtain multiple aaseqs of corresponding geneids from dictionary 1E
    geneids = []
    geneid_aaseqs = []
    for geneid in uniprotid2geneids[uniprotid].keys():
        if geneid in geneid2aaseq:
            geneids.append(geneid)
            geneid_aaseqs.append(geneid2aaseq[geneid])
    if len(geneids) == 0:
        print('No sequence available for ' + uniprotid, file=sys.stderr, flush=True)
        continue
    #    5.2 Load aaseq from the uniprodid's AlphaFold pdb file
    af_dir = '/home/bbian/Data_all/raw_data/Alpha-fold/' + AlphafoldPH
    af_pdb = af_dir + '/AF-' + uniprotid + '-F1-model_v2.pdb'
    if not os.path.isfile(af_pdb):
        print(af_pdb + ' does not exist', file=sys.stderr, flush=True)
        continue
    records = list(SeqIO.parse(af_pdb, 'pdb-seqres'))
    if len(records) > 1:
        sys.exit(af_pdb + ' contains multiple sequences')
    af_aaseq = str(records[0].seq)
    #    5.3 Compare each of aaseqs from 5.1 and aaseq from 5.2,
    #    5.4 Select the most similar aaseq from 5.1
    print('Finding the best match sequence for ' + uniprotid, file=sys.stderr, flush=True)
    best_geneid = ''
    best_identity = 0
    for i in range(len(geneids)):
        score = pairwise2.align.globalxx(geneid_aaseqs[i], af_aaseq, score_only=True)
        identity = score / len(af_aaseq)
        print('Found the sequence ' + geneids[i] + ' with the identity ' + str(identity), file=sys.stderr, flush=True)
        if identity > best_identity:
            best_geneid = geneids[i]
            best_identity = identity
            uniprotid2best_geneid[uniprotid] = best_geneid
    #    5.5 make a dictionary from uniprodid to the selected aaseq
    uniprotid2aaseq[uniprotid] = geneid2aaseq[best_geneid]
    print('Selected sequence for ' + uniprotid + ' ' + uniprotid2aaseq[uniprotid], file=sys.stderr, flush=True)

# Then, output the results with your favorite format.
# For example, for each uniprotid, output followings things as one line
for uniprotid in uniprotid2geneids.keys():
    # - uniprotid
    # - aaseq from dictionary 5.5 (not from pdb file)
    if uniprotid not in uniprotid2aaseq:
        continue
    aaseq = uniprotid2aaseq[uniprotid]
    aaseq1 = [str(i) for i in aaseq]
    aaseq1 = ','.join(aaseq1)
    # - multiple geneids (comma separated) from dictionary 1A (not from dictionary 4)
    geneids = ','.join(aaseq2geneids[aaseq])
    # - codon wave counts from dictionary 2
    if aaseq not in aaseq2codon_wave_counts:
        continue
    codon_wave_counts = aaseq2codon_wave_counts[aaseq]
    codon_wave_codons = geneid2codon_wave_codons[uniprotid2best_geneid[uniprotid]]
    codon_wave_codons = [str(i) for i in codon_wave_codons]
    codon_wave_codons = ','.join(codon_wave_codons)
    codon_wave_counts = [str(i) for i in codon_wave_counts]
    codon_wave_counts = ','.join(codon_wave_counts)
    output_values = [uniprotid,
                     aaseq1,
                     codon_wave_counts,
                     codon_wave_codons];
    print('\t'.join(output_values))

