#!/usr/bin/env python

### Goal is to count the length of each sequence and determine if divisible by 3
### Divisible by 3 == ORF in frame
### Simple remedy of trimming the 3' and then 5' (alternating) until divisible by 3
### Not intended to be the final preprocessing due to removal of sequence information
### May be accepted because core ORF still retained

from Bio import SeqIO
import os, sys
import argparse

def create_parser():
	parser = argparse.ArgumentParser(description='Count length of each nucleotide sequence to determine if overall size of ORF sequences are correct i.e., are divisible by 3; count and correction options provided.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file.")
	parser.add_argument('-o', '--output', dest="output_prefix", default="corrected", help="Output file prefix. Will be appended to the beginning of input FASTA file (i.e., 'corrected_'input.fasta)")
	parser.add_argument('--count', dest="count", action="store_true", help="Only output sequences that do not have correct frame i.e., are not divisible by 3. No output will be generated!")
	parser.add_argument('--keep-intermediate', dest="intermed_keep", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	parser.add_argument('-f', '--force', dest='force_output', action='store_true', help='By default, no output is generated if no trimming is done. This flag will force output to be created otherwise, producing a copy of the input file')
	return parser

def trim_seq(seq):
	"""
	For a given string representing a sequence, remove last character at the 3', 5' and repeat until sequence is divisible by three; thus, forcibly correcting the frame
	"""
	while len(seq) % 3 != 0:
		seq = seq[:-1] # trim last character
		if len(seq) % 3 == 0:
			return seq
		else:
			pass
		seq = seq[1:] # trim first character
		if len(seq) % 3 == 0:
			return seq
		else:
			continue

def stdout_count(failed, all):
	"""
	For a given dictionary(ies) of header: count --> send to stdout
	"""
	per_failed = (len(failed)/len(all))*100
	print(f"Original sequence lengths that require correction:"+"\n")
	for seq in failed:
		print(f"{seq}: {all[seq]}")
	print("\n"+f"{per_failed:.2f}% ({len(failed)}/{len(all)}) of sequences failed to be within frame!")

def output_seq(seqs, filename):
	"""
	For a given dictionary of header: sequence --> output to file
	"""
	with open(filename, 'w+') as out:
		for seq in seqs:
			out.write(">%s\n%s\n" %(seq, seqs[seq]))

def remove_tmp(tmp):
	if len(tmp) > 0:
		tmp = [os.remove(file) for file in tmp_files]
	return tmp

def count_frame():
	### Parser
	parser = create_parser()
	args = parser.parse_args()
	
	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	output = str(args.output_prefix)+"_"+fasta_name

	### Count each sequence, storing the length of each sequence per header]
	sequences = {}
	seq_count = {}
	failed = {} # store only failed headers, used for --count
	edit_name = 2 # ex. header2, header3, ...
	for seq_record in SeqIO.parse(fasta_file, 'fasta'):
		header = seq_record.description
		seq = str(seq_record.seq).upper()
		seq_length = len(seq)
		
		if (header in sequences) or (header in seq_count):
			header = f"{header}{edit_name}"
		sequences[header] = seq
		seq_count[header] = seq_length
		
	assert len(sequences) == len(seq_count)
	
	### Process each sequence, splicing them as needed
	for h1,h2 in zip(sequences, seq_count):
		assert h1 == h2
		if int(seq_count[h2]) % 3 != 0:
			failed[h2] = seq_count[h2]
			if not args.count: # if we just want the count, then we don't bother with trimming
				sequences[h1] = trim_seq(sequences[h1])
				assert len(sequences[h1]) % 3 == 0
				
	### Generate stdout report 
	if len(failed) == 0: # no sequences out of frame
		print("No sequences failed!")
	else:
		stdout_count(failed, seq_count)
	
	### If applicable (True, by default), generate output file with corrected sequences
	if not args.count:
		if len(failed) == 0: # no sequences out of frame
			if not args.force_output:
				print("No sequences corrected! No output will be generated!")
			else:
				output_seq(sequences, output)
				print(f"Output file: {output} generated despite no changes to sequences!")
		else:
			output_seq(sequences, output)
			print(f"Output file: {output} generated with corrected sequences!")
	
	if not args.intermed_keep:
		remove_tmp(tmp_files)

if __name__ == '__main__':
	try:
		script_path = sys.argv[0].rsplit("/",1)[0]
	except IndexError:
		script_path = os.getcwd() # script in current working directory
	start_dir = os.getcwd()
	tmp_files = []
	count_frame()
	exit()
