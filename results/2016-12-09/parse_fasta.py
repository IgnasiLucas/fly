from Bio import SeqIO
import sys

with open(sys.argv[1], "rU") as ref:
   sequences = SeqIO.parse(ref, "fasta")
   SeqIO.write(sequences, sys.stdout, "fasta")
