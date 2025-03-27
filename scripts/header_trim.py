#!/bin/env python

### Goal is to trim header and output FASTA file with headers modified
### Header is modified to be AccessionID_<bait_range> i.e. A10048.01_21-100

from Bio import SeqIO
import os, sys
from multiprocessing import Pool, Process, Lock, Queue
from itertools import repeat

try:
	input = sys.argv[1]
	threads = int(sys.argv[2])
	baits = eval(sys.argv[3])
except (IndexError, NameError):
	try:
		input = sys.argv[1]
		threads = 1
		print("Assuming input FASTA file are baits")
		baits = True
	except IndexError:
		print("""
Goal is to trim header and output FASTA file with headers modified such that it is AccessionID_<bait_range> i.e. A10048.01_21-100.

Arguments are as followed:
	header_trim.py <input_file> <num_threads (default = 1)> <baits (T/F; default =  True)>
	
	Baits argument defines whether input FASTA file contains baits generated using BaitsTools.
	""")
		exit()

def batch_iterator(iterator, batch_size): # create iterator
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

def count_seq(input):
	num = len([1 for line in open(input) if line.startswith(">")])
	return num

def split_fasta(fasta_file, num, threads):
	batch_all = [] # store split input
	record_iter = SeqIO.parse(fasta_file,"fasta")
	if threads <= 0:
		raise ValueError
	for i, batch in enumerate(batch_iterator(record_iter, int(num//threads + (num % threads > 0)))):
		i = i + 1
		if i <= threads:
			batch_all.append(batch)
		else:
			break
	return batch_all


def trim_header(batch, baits):
	header = []
	seq = []
	for seq_record in batch:
		if baits:
			header.append(str(str(seq_record.description).split(" ")[0] + "_" + str(seq_record.description).rsplit("_",1)[1]))
			seq.append(str(seq_record.seq))
		else:
			header.append(str(seq_record.id))
			seq.append(str(seq_record.seq))
	
	return header,seq
	
	#if path.exists("trimmed_" + input) is True: # do not overwrite file, concatenate
	#	with open("trimmed_" + input, "w") as out:
	#		for h,s in zip(header, seq):
	#			out.write(">" + str(h) + "\n" + str(s) + "\n")
	#else: # create file and write
	#	with open("trimmed_" + input, "w+") as out:
	#		for h,s in zip(header, seq):
	#			out.write(">" + str(h) + "\n" + str(s) + "\n")
		

def run():
	num_seq = count_seq(input)
	batch = split_fasta(input, num_seq, threads)
	if baits: print("Input FASTA file are baits")
	else: print("Input FASTA is a normal FASTA file")
	print("Split into " + str(len(batch)) + " subunit(s).")
	p = Pool(threads)
	with open(input.rsplit(".")[0] + "_trimmed." + input.rsplit(".")[1], "w") as out:
		for i,j in p.starmap(trim_header, zip(batch, repeat(baits))):
			for h,s in zip(i,j):
				out.write(">" + str(h) + "\n" + str(s) + "\n")
	#p.imap_unordered(trim_header, batch)
	#trim_header(input)

if __name__ == '__main__':
	run()
	print("Trimmed headers.")
	exit()