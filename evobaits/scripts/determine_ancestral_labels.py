#!/usr/bin/env python3

### This script is to take two FASTAs containing the same sequences, but differing headers, as input
### In addition, phylogenetic trees and JSONs produced by HyPhy Ancestral Sequences program will be taken as input
### The goal is to use the two FASTAs to create a dict object that can translate the headers from input 1 into input 2 labels
### Subsequently, we can then replace the labels in the HyPhy output to produce clearer output

### This script is an alternative approach to outright re-running the tools which may spark additional variation
### Here we do not look or tamper with sequences or even the values of the output, only the tagged label attached to the important information

import os, sys
import argparse
import re
from Bio import SeqIO # should compare sequences per header to ensure valid translation

def generate_hash(old_file, new_file):
	old = {} # seq: header
	new = {} # seq: header
	old2new = {} # old: new
	count = 1
	
	# parse individual FASTAs to extract sequences and headers
	for seq_record in SeqIO.parse(old_file, 'fasta'):
		old[str(seq_record.seq).upper()] = str(seq_record.description)
	
	if len(set([old[i] for i in old])) == 1: # all matching headers in the old, means likely a number was appended in ancestral reconstruction
		matching = True
	else:
		matching = False
	
	for seq_record in SeqIO.parse(new_file, 'fasta'):
		new[str(seq_record.seq).upper()] = str(seq_record.description)
		
	assert len(old) == len(new) # we should have equal number of sequences here
	
	# determine matching seqs and store the headers
	for seq in old:
		try:
			if matching:
				old2new[f"{old[seq]}{count}"] = new[seq]
				count+=1
			else:
				old2new[old[seq]] = new[seq]
		except IndexError: # if somehow a sequence does not exist between them
			if matching:
				old2new[f"{old[seq]}{count}"] = f"{old[seq]}{count}" # keep same header between old and new
				count+=1
			else:
				old2new[old[seq]] = old[seq]
			
	return old2new

def update_tree(file, out, headers, overwrite=False):
	intree = open(file)
	tree = ''.join(intree.readlines())
		
	for old_header in headers:
		tree = re.sub(f"{old_header}:", f"{headers[old_header]}:", tree)
		#tree.replace(f"{old_header}:", f"{headers[old_header]}:")
		
	if os.path.exists(out) and overwrite:
		os.remove(out)
	elif os.path.exists(out) and not overwrite:
		print(f"Output file: {out} already exists! Please rename or use '--overwrite'!")
		sys.exit(1)
	
	with open(out, 'w+') as outtree:
		outtree.write(tree)
	
	intree.close()
	
def update_json(file, out, overwrite=False):
	"""
	May not be needed because the tree is already a subset of the JSON.
	Ancestral nodes remain unchanged also...
	"""
	pass

def run(args):
	# Check input is valid
	if not os.path.exists(args.input_fasta_orig):
		print("Input FASTA 1 (i.e., with original unclear labels) is invalid!")
		sys.exit(1)
	if not os.path.exists(args.input_fasta_updated):
		print("Input FASTA 2 (i.e., with updated labels) is invalid!")
		sys.exit(1)
		
	old2new = generate_hash(args.input_fasta_orig, args.input_fasta_updated)
	
	
	# determine general output file name
	if args.out_pre is None:
		filename = os.path.basename(args.input_fasta_updated)
		out_pre = f"corrected_{filename}"
	else:
		out_pre = args.out_pre
	
	# determine output file prefix prior to any analysis
	if args.input_json is not None and os.path.exists(args.input_json):
		print("JSON input not supported yet! Skipping...")
		#update_json(args.input_json, args.out_pre, args.overwrite)

	if args.input_tree is not None and os.path.exists(args.input_tree):
		out_path = os.path.join(args.out_dir, f"{out_pre}.tree")
		update_tree(args.input_tree, out_path, old2new, args.overwrite)

def create_parser():
	parser = argparse.ArgumentParser(description='Translation script to update labels in HyPhy AncestralSequences JSON and TREE output as generated via ViralSTEP clustering.')
	parser.add_argument('-i1', '--input1', dest="input_fasta_orig", required=True, help="Input FASTA file with unclear headers. This should be the FASTA from the early ViralSTEP clustering.")
	parser.add_argument('-i2', '--input2', dest="input_fasta_updated", required=True, help="Input FASTA file with updated or taxonomically identified headers. This should be the output from 'determine_taxonomy.py'.")
	parser.add_argument('-it', '--input-tree', dest='input_tree', help="Input tree file to update labels.")
	parser.add_argument('-ij', '--input-json', dest='input_json', help="NOT SUPPORTED YET! Input JSON file to update labels.")
	parser.add_argument('-od', '--output-dir', dest="out_dir", default=os.getcwd(), help="Output directory for corrected files. Default in the current working directory.") 
	parser.add_argument('-op', '--output-prefix', dest='out_pre', help="Output prefix for resulting translated files. Default will retain input 2 filename with appended 'corrected'.")
	#parser.add_argument('-f', '--field', dest='field', type=int, default=1, help='Field position with "_" delimiter where ORF name is to be expected if filename is altered.')
	parser.add_argument('--overwrite', action='store_true', help='If output file(s) already exist, overwrite them.')
	
	return parser
	
if __name__ == '__main__':
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	tmp_files = []
	
	try:
		script_dir = script_path.rsplit("/", 1)[0]
	except IndexError:
		script_dir = os.getcwd() # script in current working directory
		
	parser = create_parser()
	args = parser.parse_args()
	
	run(args)
	print("Correction completed successfully!")
	sys.exit(0)