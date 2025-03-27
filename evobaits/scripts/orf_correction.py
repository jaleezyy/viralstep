#!/usr/bin/env python

### Goal is to mitigate as many errors with FASTA format
### Errors include odd characters and stop codon characters that error out in HyPhy
### Will act as a preprocessing step prior to alignment

from Bio import SeqIO
import os, sys
import argparse
from convert import detect_content # convert.py

def create_parser():
	parser = argparse.ArgumentParser(description='Given nucleotide or protein ORF sequences, examine for consistent characters corresponding to the content and remove stop codon traces, if desired.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file. Will auto-detect content.")
	parser.add_argument('-o', '--output', dest="output_prefix", default="corrected", help="Output file prefix. Will be appended to the beginning of input FASTA file (i.e., 'corrected_'input.fasta)")
	parser.add_argument('-c', '--content', dest="input_content", choices=['nucl', 'prot'], default=None, help="Sequence content for input sequences. Choose between 'nucl' or 'prot'. Leave blank for auto-detection of content per sequence.")
	parser.add_argument('--offset', dest="offset", default=9, help="Determine how many characters/codons will be examined per sequence at the 5' and 3' ends to validate start/stop codons. Default = 9")
	parser.add_argument('--all-stops', dest='remove_all_stops', action='store_true', help="Remove all traces of stop codons")
	parser.add_argument('--keep-intermediate', dest="intermed_keep", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	return parser

			
def check_seq(seq, valid_subset, all_valid, remove=False):
	"""
	Given a sequence (partial or otherwise) and a list of valid characters, check and trim as needed, returning the corrected sequences.
	
	If start or stop, we don't want to trim characters unless they generally are invalid (i.e., are not a valid key). Not every codon or character will match the subset of keys denoted for their respective section (start or stop)
	
	The 'type' variable will be used to determine the split pattern (1 for AA; 3 for nucleotide aka codon sequences)
	
	Portion is an internal variable highlighting whether we're examining the 5', 3' ends or the remaining sequence (default)
	"""
	# check content to determine splitting of provided sequence
	# alternative: use valid (length of value(s)) to determine this split
	if len(valid_subset[0]) == 1: # amino acid 
		content = 1
	elif len(valid_subset[0]) == 3: # nucleotide codon
		content = 3
	else:
		content = 1
		
	# split sequence by specified chunk
	sequence = [seq[i:i+content] for i in range(0, len(seq), content)]

	# for each chunk, compare with valid keys to ensure nothing is amiss
	for chunk in sequence:
		if chunk in valid_subset:
			# if chunk matches, but we want to remove it (i.e., stop codon), then we do it here
			if remove:
				index = sequence.index(chunk)
				sequence[index] = "-"*content
			else:
				continue
		else: 
			# determine if chunk contains a valid character or codon, ignore if sort
			# only replace with gaps if an error occurs checking against all potetnially valid keys
			if str(chunk).upper() in all_valid:
				continue
			elif (content == 3) and (any(nucl in str(chunk).upper() for nucl in ['A', 'C', 'T', 'U', 'G'])):
				continue
			else:
				index = sequence.index(chunk)
				sequence[index] = "-"*content
	
	final_seq = ''.join(map(str,sequence))
	
	return final_seq

def trim_seq(seq, offset=9, type=None, remove_all=False):
	"""
	For a given string representing a sequence, obtain the first and last <offset> (AA or <offset>*3 nucleotides, subject to change) to correct start and stop portions (i.e, trim stop codons, validate methionine, etc). From there, check for valid characters throughout the rest, trim as needed.
	
	
	Return corrected sequence
	"""
	### TODO: Use Biopython to generate this table 
	### Bacterial, archaeal and plant plastid code (including prokaryotic viruses)
	tranl_table_11={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", 
        "TGT": "C", "TGC": "C", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
	stop_codons=["TAA", "TAG", "TGA"]
	start_codons=["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"] # alternative start codons
	
	# check length of sequence to ensure offset won't error
	if len(seq) < offset:
		offset = len(seq)
	
	# determine sequence content
	if type is None:
		type = detect_content(seq)
	if type == "nucl":
		offset = offset*3
	#print(f"{type} detected!")
	first = seq[:offset] # 5' only
	last = seq[-offset:] # 3' only
	remaining = seq[offset:-offset] # trimmed 5' and 3'
	
	# separate sequence into three parts: 5'-end, 3'-end, and remaining sequence
	if type == "nucl": # keep nucleotide comparisons
		valid_keys = [codon for codon in tranl_table_11]
		valid_stop = stop_codons
		valid_start = start_codons
	elif type == "prot": # convert to amino acid keys
		valid_keys = list(tranl_table_11.values())
		valid_stop = ["*"] # defaut stop codon character
		valid_start = [tranl_table_11[s] for s in start_codons]
	
	# check start codon
	first = check_seq(first, valid_start, valid_keys)
	if remove_all:
		first = check_seq(first, valid_stop, valid_keys, True)
	
	# check stop codons
	last = check_seq(last, valid_stop, valid_keys, remove_all)
	
	# check remaining sequence
	remaining = check_seq(remaining, valid_keys, valid_keys)
	if remove_all:
		remaining = check_seq(remaining, valid_stop, valid_keys, True)
	
	# merge resulting sequence
	final_seq = str(first+remaining+last)
	
	return final_seq

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

def correct_frame():
	### Parser
	parser = create_parser()
	args = parser.parse_args()
	
	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	output = str(args.output_prefix)+"_"+fasta_name
	offset = int(args.offset) 
	edit_name = 2 # use to modify matching header, increment by 1
	
	### Collect each sequence
	sequences = {}
	seq_count = {}
	for seq_record in SeqIO.parse(fasta_file, 'fasta'):
		header = seq_record.description
		seq = str(seq_record.seq).upper()
		
	### Process each sequence, splicing them as needed
		seq = trim_seq(seq, offset, args.input_content, args.remove_all_stops)
	
	### Store final header-sequence pair for output
		if header in sequences:
			header = f"{header}{edit_name}" # ex. header2, header3, ...
			edit_name+=1
		sequences[header] = seq
	
	### If applicable (True, by default), generate output file with corrected sequences
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
	correct_frame()
	exit()
