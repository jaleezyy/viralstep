#!/bin/env python

from Bio import SeqIO
import os, argparse, sys, itertools 
#import subprocess

# Script aims to take a single FASTA file and split it between a number of files
# Total number of output files depends on the desired number of sequences per file
# Example: split a FASTA containing 10 sequences by 1 --> will generate 10 FASTA files with 1 sequence each
# Example 2: split a FASTA containing 100 sequences by 5 --> will generate 20 FASTA files with 5 sequences each

#userParameters = sys.argv[1] # input fasta file that will be filtered
#split = sys.argv[2]
#final = userParameters+"_split"

def parser():
	args = []
	try:
		input_fasta = sys.argv[1]
		if os.path.isfile(input_fasta):
			args.append(input_fasta)
		else:
			raise IndexError
		batch_size = int(sys.argv[2])
		args.append(batch_size)
		return args
	except IndexError:
		return args


def batch_iterator(iterator, batch_size):
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = next(iterator)
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch


def split_fasta(fasta_file, split):
	#os.system("mkdir %s" % (final))
	#os.system("cd ./%s'_split'" % (userParameters))
	record_iter = SeqIO.parse(fasta_file,"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, int(split))):
		filename = "group_0%i_"  % (i + 1) + os.path.basename(fasta_file)
		with open(filename, "w") as handle:
			count = SeqIO.write(batch, handle, "fasta")
		print("Wrote %i records to %s" % (count, filename))
	
#split_fasta(userParameters, split)
#os.system("cd -")
#os.system("mv group_*_* ./%s" %(final))
	
#print("Done.")
if __name__ == '__main__':
	HELP = """
Script aims to take a FASTA file as input and create split FASTA output with a specified number of sequences per file as a maximum.

Usage:
        split_fasta.py <input_fasta> <batch_size>
"""
	args = parser()
	print(args)
	if len(args) == 0 or len(args) > 2:
		print(HELP)
	else:
		split_fasta(args[0], args[1])
