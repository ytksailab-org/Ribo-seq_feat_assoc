import sys
import csv
import re
from Bio import SeqIO
from Bio.Seq import Seq

#usage
#python nuc2codon.transcriptome.mapping.py yeast.transcriptome.fa yeast_Riboseq_calbrited.nucle.wave yeast_Riboseq_calbrited.codon.wave
input_fasta = sys.argv[1]  # yeast.transcriptome.fa
input_nuc_wave = sys.argv[2]  #yeast_Riboseq_calbrited.nucle.wave
out_codon_wave = sys.argv[3]  #yeast_Riboseq_calbrited.codon.wave

csv.field_size_limit(100000000000)
record_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
#print(str(len(record_dict)) + ' sequences loaded')
#print(str(record_dict))

for name in record_dict.keys():
    record = record_dict[name]
    #print(name + ' has the sequence length ' + str(len(record)))

with open(out_codon_wave, mode='w') as out_f:
    writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')

    with open(input_nuc_wave) as f:
        # Assume each line has tab-separated values.
        reader = csv.reader(f, delimiter='\t')
        for row in reader:  # row here is list
            genename = row[0]
            #print(genename)
            if genename not in record_dict:  # ensure the name is same
                sys.exit(genename + ' does not exist in fasta')

            cds_seq_old = record_dict[genename].seq
            descriptons = record_dict[genename].description
            if "CDS" in descriptons:
                m= descriptons.split("CDS=")[1]
                start= m.split('-')[0]
                stop = m.split('-')[1]
            else:
                continue
            cds_seq = cds_seq_old[int(start) - 1: int(stop)]
            #cds_seq = str(cds_seq)
            if len(cds_seq) % 3 == 1:
                cds_seq = cds_seq + Seq('NN')
            elif len(cds_seq) % 3 == 2:
                cds_seq = cds_seq + Seq('N')

            protein_seq = cds_seq.translate()
            print('Eextracted nucleotide sequence:')
            print(cds_seq)
            print('Translated amino acid sequence:')
            print(protein_seq)
            #Report if translated CDS does not end with stop codon *
            # or contains internal stop codons.
            # This can happen for pseudogenes.
            if not protein_seq[-1] == '*':  # one char from the end
                print('info: translated CDS does not end with *')
            if '*' in protein_seq[:-1]:  # exclude the last character
                print('info: translated CDS contains internal stop codons *')
            if protein_seq[-1] == '*' and '*' not in protein_seq[:-1]:

                nuc_wave_dict = {}  # create empty dict
                nuc_wave_positions = row[1].split(',')  # split the str into the list using ','
                nuc_wave_counts = row[2].split(',')
                if not len(nuc_wave_positions) == len(nuc_wave_counts):  # should be the same, validate check
                    sys.exit('wrong length of nucleotide-wise wave')
                for i in range(len(nuc_wave_positions)):  # for each index in the list
                    nuc_wave_dict[int(nuc_wave_positions[i])] = int(nuc_wave_counts[i])  # creat the dict key= value, non-countinous

                # Convert nucleotide-wise wave into codon-wise wave
                codon_wave_positions = []  # creat 4 empty list
                codon_wave_counts = []
                codon_wave_codons = []
                codon_wave_aas = []
                # If plus strand, read the wave in 5-->3 direction
                # If minus strand, read the wave in 3-->5 direction

                wave_start = int(start)-1
                wave_stop = int(stop)   # range(1:10)
                wave_dir = +1  # direction


                # For each codon ...
                idx = 0  # number of codon processed
                for i in range(wave_start, wave_stop, wave_dir * 3):  # stipe
                    p1 = i  # p  position in the codon
                    p2 = i + wave_dir
                    p3 = i + wave_dir + wave_dir
                    c1 = 0  # set default value of empty reads position
                    c2 = 0
                    c3 = 0
                    if p1 in nuc_wave_dict:
                        c1 = nuc_wave_dict[p1]
                    if p2 in nuc_wave_dict:
                        c2 = nuc_wave_dict[p2]
                    if p3 in nuc_wave_dict:
                        c3 = nuc_wave_dict[p3]
                    codon_wave_positions.append(str(p1))
                    codon_wave_counts.append(str(c1 + c2 + c3))
                    codon_wave_codons.append(str(cds_seq[idx * 3:idx * 3 + 3]))  # three character codons sequence
                    codon_wave_aas.append(str(protein_seq[idx]))  # one amino acid character in this position
                    idx += 1
                # Length of codon wave and protein sequence should be the same.
                if not len(codon_wave_positions) == len(protein_seq):
                    sys.exit('inconsistent sequence length')
                # Format output line.

                output_values = ["NC_000001.11",
                                 genename,
                                 ','.join(codon_wave_positions),
                                 ','.join(codon_wave_counts),
                                 ','.join(codon_wave_codons),
                                 ','.join(codon_wave_aas)];
                writer.writerow(output_values)  # print('\t'.join(output_values))  # use tab to combine the list


