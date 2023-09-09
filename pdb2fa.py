import sys
from Bio import SeqIO
af_pdb = sys.argv[1]
seqs = sys.argv[2]

records = list(SeqIO.parse(af_pdb, 'pdb-seqres'))
#print(records)
SeqIO.write(records[0], seqs, "fasta")
