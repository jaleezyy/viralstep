#!/bin/env python

### Goal is to pull random sequences from a given FASTA file
### 4 arguments - input file, db (accID, taxID, baltimore, virus family) file, output file, length of pull (how many elements?) 

import sys, os
import argparse
from Bio import SeqIO
from random import sample
from collections import defaultdict

# try:
	# if "help" in [x.lower() for x in sys.argv]: raise IndexError
	# input = sys.argv[1]
	# db_file = sys.argv[2]
	# output = sys.argv[3]
	# pull_length = sys.argv[4]
# except IndexError:
	# print("""
# Arguments are as followed:
	# random_seq_pull.py <input_file> <db_file> <output_file> <length of randomized pull>
	
# Note: db_file is a CSV file with the following columns: accession ID, taxonomic ID, balitmore classifier, virus family
	# """)
	# exit()
	
def create_parser():
	parser = argparse.ArgumentParser(prog='viralstep normalize', description='Generate a subset from a given bait set with a set maximum on number of sequneces per virus family.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file")
	parser.add_argument('-db', '--database', dest="db_file", required=True, help="Input database file in CSV format")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=False, default="subset.fasta", help="Output FASTA file. Default = 'subset.fasta'")
	parser.add_argument('-s', '--size', dest="pull_length", required=True, help="Maximum number of sequences per virus family")
	parser.add_argument('--generatedb', dest="subset_db", action="store_true", help="Generate CSV DB file based on subset FASTA")
	return parser

def load_reference(db):
	print("Reading reference database...")
	syndrome_list = set() # list of individual virus families
	accID2taxID = defaultdict(list) # taxID: [accIDs]
	taxID2family = {} # normal dictionary object since one taxID = one family
	family2baltimore = {} # family: baltimore classification (i.e. dsDNA, +ssRNA, etc.)

	with open(db) as f:
		for line in f:
			syndrome_list.add(line.split(",")[-1])
			accID2taxID[line.split(",")[1]].append(line.split(",")[0])
			taxID2family[line.split(",")[1]] = line.split(",")[-1]
			family2baltimore[line.strip("\n").split(",")[-1]] = line.split(",")[-2]
	
	return [syndrome_list], accID2taxID, taxID2family, family2baltimore
	#return syndromes, accID2taxID, taxID2family, family2baltimore

def random_sample(seqs, pull):
	try:
		randomized = sample(seqs, int(pull))
	except ValueError: # number of seqs < pull
		return seqs
	return randomized

def run(input, syndromes, accID2taxID, taxID2family, pull_length):
	#syndromes, accID2taxID, taxID2family = load_reference(db_file)
	seq_per_syn = defaultdict(list)

	### Group individual FASTAs by virus family
	for seq_record in SeqIO.parse(input, 'fasta'):
		print("Identifying virus family for %s" %(seq_record.id))
		for x,y in accID2taxID.items():
			# isolates AccessionID from any sequence header (including post BaitsTools)
			# Interpreting bait header
			if seq_record.description.split(" ")[0].rsplit("_",1)[0] in y: 
				virus_family = str(taxID2family[x]).split("\n")[0] # virus family
			# Interpreting regular accession ID
			# Opposite true if db file contains bait headers
			elif seq_record.description.split(" ")[0] in y: 
				virus_family = str(taxID2family[x]).split("\n")[0] # virus family
			else:
				continue
		#virus_family = taxID2family[search_term]
		print(str(seq_record.id)+": "+str(virus_family))
		seq_per_syn[virus_family].append(seq_record)

	### Randomly pull from each virus family (list from list)
	print("Pulling random sequences per virus family")
	subsample = {}
	for fam in seq_per_syn:
		if "unclassified" in str(fam).lower(): continue
		if len(seq_per_syn[fam]) > 0:
			subsample[fam] = random_sample(seq_per_syn[fam], pull_length)

# with open(input, 'r') as seq_record:
	# seqs = SeqIO.parse(seq_record, 'fasta')
	# # parse list of randomized samples and re-run if needed
	# seq_pull = random_sample(input)
	# seq_syndrome = taxID2family[accID2taxID[seq_pull.id]]
	# try:
		# if len(seq_per_syn[seq_pull.id]) < pull_length:
			# seq_per_syn[seq_syndrome].append(seq_pull)
	# except KeyError:
		# seq_per_syn[seq_syndrome].append(seq_pull)
	
	# if not all(len(seq_per_syn[i]) == pull_length for i in syn_iter):
		# continue
	# else:
		# break
	
	return subsample
	
### OUTPUT
def write_output(output, subsample):
	print("Outputting select sequences")
	with open(output, 'w') as out:
		for seq in subsample:
			for seq_record in subsample[seq]:
				out.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")

def update_db(subset_fasta, accID2taxID, taxID2family, family2baltimore):
	output_path = os.path.dirname(subset_fasta)
	output_file_name = str(os.path.basename(subset_fasta)).rsplit(".", 1)[0]
	final_db = os.path.join(output_path, "db_%s.csv" %(output_file_name))
	print("Generating DB file: %s" %(final_db))
	
	with open(final_db, 'w+') as db_out:
		for seq_record in SeqIO.parse(subset_fasta, 'fasta'):
			for x,y in accID2taxID.items():
				# isolates AccessionID from any sequence header (including post BaitsTools)
				# Interpreting bait header
				if seq_record.description.split(" ")[0].rsplit("_",1)[0] in y: 
					accID = seq_record.description.split(" ")[0].rsplit("_",1)[0] ### convert to seq_record.id
					taxID = x
					virus_family = str(taxID2family[x]).split("\n")[0] # virus family
					balt_class = str(family2baltimore[virus_family])
				# Interpreting regular accession ID
				# Opposite true if db file contains bait headers
				elif seq_record.description.split(" ")[0] in y:
					accID = seq_record.description.split(" ")[0] ### convert to seq_record.id
					taxID = x
					virus_family = str(taxID2family[x]).split("\n")[0] # virus family
					balt_class = str(family2baltimore[virus_family])
				else:
					continue
			db_out.write("%s,%s,%s,%s\n" %(accID, taxID, balt_class, virus_family))

	pass

if __name__ == '__main__':
	parser = create_parser()
	args = parser.parse_args()
	syndromes, accID2taxID, taxID2family, family2baltimore = load_reference(args.db_file)
	subset = run(args.input_fasta_file, syndromes, accID2taxID, taxID2family, args.pull_length)
	write_output(args.output_fasta_file, subset)
	if args.subset_db: 
		update_db(args.output_fasta_file, accID2taxID, taxID2family, family2baltimore)
	else: 
		pass
	print("Pull complete!")
	exit()