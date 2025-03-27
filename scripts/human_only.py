#!/bin/env python

### Goal is to create a script to subset an existing bait set to only those baits that target human-infecting viruses
### Virus family independent (i.e., respiratory syndromic baits filtered for human-host respiratory viruses)
### Uses virustohost DB (run wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv to download)

import os, sys, subprocess, shutil
import argparse
import datetime, time
from Bio import SeqIO
import urllib.request


def create_parser():
	parser = argparse.ArgumentParser(prog="viralstep humanwgc", description='Filter given existing baitset to a subset targetting human-host infecting viruses.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file. Should conform to the header structure AccID_##-##.")
	parser.add_argument('-db', '--database', dest="virus2hostdb", required=False, help="Virus-Host DB file from ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv (ViralSTEP default: will attempt to search within VHost-Classifier).")
	parser.add_argument('-l', '--list', dest="taxid_list", required=True, help="Input CSV file containing individual bait accessionID, taxID, baltimore classification, virus family. At minimum, the accessionID and taxID columns are required.")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="Output fasta file")
	parser.add_argument('--keep-intermediate', dest="output_pre", action="store_true", required=False, help="Keep intermediate temporary files.") ### will delete at end of program otherwise
	return parser

def count_seq(input):
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq

def download_db(dir):
	out_db = os.path.join(dir, "virushostdb.tsv")
	if os.path.exists(out_db):
		print("Existing Virus-Host DB file found! Updating...")
	print("Downloading Virus-Host DB file to %s/virushostdb.tsv" %(os.path.abspath(dir)))
	with urllib.request.urlopen("ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv") as response, open(out_db, 'wb') as out_file:
		shutil.copyfileobj(response, out_file)
	#tmp_files.append("virustohost.tsv")
	return out_db

### TODO: Add module to generate accessionID, taxID from given FASTA file if no list given 
### Paired with validate_taxid functon below
def generate_taxID(input_fasta, db):
	pass

### TODO: Validation of TaxID list can be done in two ways: length of object i.e., # of seqs and IDs match and if md5sum is consistent between generated list and given list 
### md5sum check would validate both content and order of list
### In both cases (no list and given list) a taxID list will be generated (which is why validation step has optional flags to skip and only perform validation)
### This ensures a normal run with defaults only (which assumes no given list) is most accurate (but moddable)
def validate_taxid(input_list, generated_list, num_seq):
	pass

### TODO: Recycled old function - update to dataframe interpretation
def load_references(taxid_list, db):
	#accID2taxID = defaultdict(list)
	print("Loading reference files.\n")
	accID2taxID = {}
	human_host = set()
	
	with open(taxid_list) as f:
		for line in f:
			#accID2taxID[line.split(",")[1]].append(line.split(",")[0])
			if str(line.split(",")[0]).count("_") > 1:
				accID2taxID[line.split(",")[0].rsplit("_",1)[0]] = line.split(",")[1] # may be a memory hog
			else:
				accID2taxID[line.split(",")[0]] = line.split(",")[1]
	
	with open(db) as h:
		for line in h:
			chunk = line.split('\t')
			virus_id = chunk[0]
			human_id = chunk[7]
			#human_taxid = 9606
			if human_id == "9606":
				human_host.add(virus_id)
			else:
				continue
		
	return accID2taxID, human_host

### Modify for multiprocessing
def writer(output, queue_filter, stop): # generate output through queues
	if output not in (".fa", ".fasta"): # ensure final output is a .fasta file
		output = output + ".fasta"
	else:
		pass
### Output core file (filtered list)		
	with open(output, "w+") as out:
		while True:
			line_filter = queue_filter.get()
			if stop in line_filter:
				return
			else:
				out.write(">" + str(line_filter.description) + "\n" + str(line_filter.seq).upper() + "\n")

def filter_human(input_fasta, accID2taxID, human_host, out_file):
	print("Filtering baits.\n")
	### Parameters
	input_num_seq = count_seq(input_fasta)
	kept_baits = 0
	missing_tax = 0
	missing_id = set()
	num_missing_id = 0
	fasta_name = str(os.path.basename(input_fasta))
	out_name = str(os.path.basename(out_file))
	
	with open(out_file, 'w+') as out:
		for seq_record in SeqIO.parse(input_fasta, 'fasta'):
			in_acc = seq_record.description.split(" ")[0].rsplit("_",1)[0] # TODO: add check for NC_Acc_1-80 type format
			# for tax,acc in accID2taxID.items():
				# if in_acc not in acc: continue
				# else:
					# if tax in human_host:
						# pass
			try:
				in_tax = accID2taxID[in_acc]
				if in_tax in human_host:
					print("Human host: Bait %s" %(in_acc))
					out.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq).upper() + "\n")
					kept_baits+=1
			except KeyError:
				try:
					in_acc = seq_record.description.split(" ")[0] # regular accession (i.e., NC_100000)
					in_tax = accID2taxID[in_acc]
					if in_tax in human_host:
						print("Human host: Bait %s" %(in_acc))
						out.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq).upper() + "\n")
						kept_baits+=1
				except KeyError:
					print("Missing host: Bait %s" %(seq_record.description.split(" ")[0]))
					missing_id.add(in_acc)
					missing_tax+=1
					continue
	
	### Collect Accession IDs with missing TaxIDs (i.e., go back and check TaxID list and re-run)
	with open("skipped_from_%s_to_%s.txt" %(fasta_name, out_name), 'w+') as no_out:
		for id in missing_id:
			no_out.write(str(id) + "\n")
	
	count_baits = input_num_seq - kept_baits # number of baits removed from starting bait set
	
	if len([1 for line in open("skipped_from_%s_to_%s.txt" %(fasta_name, out_name))]) == 0:
		tmp_files.append("skipped_from_%s_to_%s.txt" %(fasta_name, out_name))
	else:
		num_missing_id = len(missing_id)
	
	return count_baits, missing_tax, num_missing_id

def clean(tmp_files):
	for file in tmp_files:
		if os.path.exists(file): os.remove(file)

### TODO: Add multiprocessing component
def run():
	### Parser
	parser = create_parser()
	args = parser.parse_args()
	
	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	out_file = os.path.abspath(args.output_fasta_file)
	
	### Download virustohost DB
	if args.virus2hostdb is None:
		virus2hostdb_path = download_db(os.getcwd())
		tmp_files.append(virus2hostdb_path)
	else:
		print("Using input virushostdb file")
		virus2hostdb_path = args.virus2hostdb
	
	### Load accID --> taxID + human-host viruses
	accID2taxID, human_host = load_references(args.taxid_list, virus2hostdb_path)
	
	### Filter bait set for human-host targetting viruses
	count_baits, missing_tax, missing_acc = filter_human(fasta_file, accID2taxID, human_host, out_file)
	
	### Clean up intermediate files
	if not args.output_pre: clean(tmp_files)
	else: print("\nKeeping intermediate files.")

	return count_baits, missing_tax, missing_acc

if __name__ == '__main__':
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	tmp_files = []
	try:
		script_path = sys.argv[0].rsplit("/",1)[0]
	except IndexError:
		script_path = os.getcwd() # script in current working directory
	start_dir = os.getcwd()
	num_baits, missing_baits, missing_acc = run()
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	print("\nRemoved " + str(num_baits) + " bait(s).\nNumber of baits skipped: " + str(missing_baits) + "\nNumber of Accession IDs with missing Taxonomic IDs: " + str(missing_acc) + "\n")
	print("Start: " + start)
	print("Finish: " + finish)
