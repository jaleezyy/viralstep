#!/bin/env python

from Bio import SeqIO
import sys

# try:
	# userParameters = sys.argv[1]
	# output = userParameters.rsplit(".",1)[0]
# except Exception:
	# exit("Requires single input .FASTA file.")

def fasta2bed(input_fasta):
	output_bed = str(input_fasta).rsplit(".",1)[0]
	with open(input_fasta) as in_f, open(output_bed+'.bed','w') as out_f:
		for record in SeqIO.parse(in_f, 'fasta'):
			out_f.write('{}\t0\t{}\n'.format(record.id, len(record)))

if __name__ == '__main__':
	fasta2bed(sys.argv[1])
