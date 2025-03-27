#!/usr/bin/env python

### Goal is to have function dedicated to detection and conversion of nucleotide <--> protein
### Codon bias calculations will determine probability of reverse translation
### Note: these scripts operate with a single sequence as input

from Bio import SeqIO
from Bio.Seq import Seq
import sys, os, argparse


def create_parser():
	parser = argparse.ArgumentParser(description='Convert given FASTA file between nucleotide and protein sequences.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file. Will auto-detect for nucleotide or amino acid content")
	parser.add_argument('-c', '--codon', dest="input_codon_freq", default=None, help="Input codon frequency table. Required for protein -> nucleotide conversion")
	parser.add_argument('-o', '--output', dest="output_prefix", default="converted", help="Output file prefix. Will be appended to the beginning of input FASTA file (i.e., 'converted_'input.fasta)")
	parser.add_argument('--codon-only', dest="codon_only", action="store_true", help="Run codon frequency calculations only for given NUCLEOTIDE FASTA")
	parser.add_argument('--convert-only', dest="convert_only", action="store_true", help="Run sequence conversion only")
	parser.add_argument('--keep-intermediate', dest="intermed_keep", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	return parser

def count_seq(input):
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq

def detect_content(input_seq):
	### Assume nucleotide (nucl) by default, caluclate freq of each nucleotide
	### If nucleotide sequence, added frequency should be 1 (will assume RNA 'U')
	type = "nucl"

	sequence = str(input_seq).upper()
	num_char = len(sequence)
	a = sequence.count("A")/num_char
	c = sequence.count("C")/num_char
	g = sequence.count("G")/num_char
	t = sequence.count("T")/num_char
	u = sequence.count("U")/num_char
	# likely 0 more often than not (unless AA)

	if a+c+g+t+u != 1:
		type = "prot"
		
	return type

def codon_usage(file, out='codon_index.cai'):
	from CAI import CAI, relative_adaptiveness # only applicable at this function
	### Can only work with nucleotide FASTA file
	#output = open(out, w+)
	
	# Calculate weighted references (remove redundant calculations)
	weights = relative_adaptiveness(sequences=[str(seq_record.seq).upper() for seq_record in SeqIO.parse(file, 'fasta')])
	
	for seq_record in SeqIO.parse(file, 'fasta'):
		type = detect_content(seq_record.seq)
		if type != "nucl":
			return "Error" # error out
		else:
			index = CAI(str(seq_record.seq).upper(), weights=weights)
			print(index)
	
	#output.close()

def convert_nucl_to_prot(nucl):
	### Input should be Seq(nucl_seq) i.e., a SeqIO object
	### Output will be Seq(prot_seq)
	sequence = Seq(str(nucl).upper())
	prot = sequence.translate(table=11, stop_symbol='*', to_stop=False, cds=False, gap='-')
	#prot = sequence.translate()

	return prot
	
def convert_prot_to_nucl(prot, codon):
	### Input should be Seq(prot_seq) i.e., a SeqIO object
	### Output will be Seq(nucl_seq)
	sequence = Seq(str(prot).upper())
	
	
	pass

def convert_seq(input_seq, codon=None):
	### Take single sequence, determine starting content (i.e., nucleotide or protein)
	### Depending on type, convert to opposite
	type = detect_content(input_seq)
	
	if type == "nucl":
		converted = convert_nucl_to_prot(input_seq)
	elif (type == "prot") and (codon is not None):
		converted = convert_prot_to_nucl(input_seq, codon)
	else:
		return "Unable to determine sequence content! If providing an amino acid sequence, ensure codon bias provided!"
		
	return converted

def convert():
	# tmp_files = []
	### Parser
	parser = create_parser()
	args = parser.parse_args()

	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	output = str(args.output_prefix)+"_"+fasta_name
	start_num = count_seq(fasta_file)
	
	print("Number of sequences in %s is %s" %(fasta_name, start_num))
	
	### Take FASTA file and run codon usage calculations + conversions as needed
	# Calculate codon usage cumulatively
	if args.convert_only:
		codon_table = None
	elif args.input_codon_freq is not None:
		codon_table = args.input_codon_freq
	else:
		codon_table = codon_usage(fasta_file)
		if args.codon_only: return None
	
	# Run conversion after calculating codon usage
	with open(output, 'w+') as out:
		for seq_record in SeqIO.parse(fasta_file, 'fasta'):
			header = str(seq_record.description)
			new_seq = convert_seq(seq_record.seq, codon_table)
			if new_seq.startswith("Unable"): 
				continue
			out.write(">%s\n%s\n" %(header, new_seq)) 
	
	end_num = count_seq(output)
	
	print("Number of sequences converted from %s to %s is %s" %(fasta_name, output, end_num)) 

	# if not args.intermed_keep:
		# tmp_files = [os.remove(tmp) for tmp in tmp_files]
	
if __name__ == '__main__':
	try:
		script_path = sys.argv[0].rsplit("/",1)[0]
	except IndexError:
		script_path = os.getcwd() # script in current working directory
	start_dir = os.getcwd()
	convert()
	exit()