#!/bin/env python

### Goal is to run VHost-Classifier and subsequent filtration script (already built)
### TODO: Merge this script and parser (to actually filter sequences) script into one file
### TODO: See what I can do to make the script more flexible in terms of pathing

import os, sys, shutil
from Bio import SeqIO
import subprocess
import glob
import argparse

"""
try:
	if "help" in sys.argv: raise IndexError
	input_fasta = sys.argv[1]
	output_fasta = sys.argv[2]
	input_csv = sys.argv[3]
except IndexError:
	exit('''
Wrapper script to execute VHost-Classifier filtration of sequences. 
Currently set to filter for Mammalian virus sequences.

	vhost_filter.py <input_fasta> <output_fasta> <db_csv>

Note: db_csv refers to a comma-separated file with accessionID,taxID,baltimore,virus_family.

''')
"""
def create_parser():
	parser = argparse.ArgumentParser(prog="viralstep vhost_filter", description='Filter sequences based on infection host. Runs VHost-Classifier + subsequent filter for mammalian virus sequences.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="input fasta file")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="output fasta file")
	parser.add_argument('-c', '--classify', dest="input_db_csv", required=True, help="input .csv with accessionID, taxID, baltimore, virus_family")
	return parser

### Parse db into {accession: (taxID, virus_family)} + output taxid list
def parse_db(input_csv):
	print("Parsing provided DB file.\n")
	db = {}
	for line in open(input_csv):
		accession = line.split(",")[0]
		taxid = line.split(",")[1]
		baltimore = line.split(",")[-2]
		virus_family = line.split(",")[-1]
		db[accession] = (taxid, baltimore, virus_family)
	# this file is the only exception that DOES NOT need to be found within VHost-Classifier
	# currently hard-coding this to VHost-Classifier directory until code pathing more flexible
	with open(os.path.join(directory, "VHost-Classifier/.vhost_filter.taxid"), 'w+') as tax_out:
		for key,value in db.items():
			tax_out.write(str(value[0]) + "\n")
	return db
	
### Run VHost-Classifier
### Currently VHost-Classifier only works from within its directory
### TODO: Rewrite VHost-Classifier code such that it is more flexible in pathing
### TODO: Use os.system('locate') to find required files?
### Add optional update flag to virushostdb 
def vHost(start_dir):
	print("Running VHost-Classifier.\n")
	os.chdir(os.path.join(directory, "VHost-Classifier"))
	# VHost_Classifier.py <list_taxids> <vhost_db> <output_dir>
	subprocess.run(['python', 'VHost_Classifier.py', '.vhost_filter.taxid', 'virushostdb.tsv', 'vhost_filter']) #, stdout=subprocess.PIPE, universal_newlines = True)
	os.chdir(start_dir)
	
### Filter based on VHost-Classifier output_dir
### TODO: Consider importing function instead of executing the script (2nd choice)
### TODO: Just migrate vhost_parse.py key function and re-work variables 
def parse_vHost(input_fasta, output_fasta):
	print("\nParsing VHost-Classifier output for Mammalian viruses.\n")
	print("(1/3) Combining all mammalian virus results from VHost-Classifier\n")
	mammal = os.path.join(directory, 'VHost-Classifier/vhost_filter/Virus/Host Assigned/Eukaryote/Chordata/Mammalia')
	vhost_csv = os.path.join(directory, 'VHost-Classifier/.vhost_filter_mammal.csv') 
	with open(vhost_csv, 'a') as out_mammal:
		for filename in glob.glob(os.path.join(mammal, '*.csv')):
			if "counts.csv" in str(filename).lower():
				continue
			else:
				pass
			with open(filename) as read_file:
				for line in read_file:
					out_mammal.write(line)
	print("(2/3) Outputting filtered sequences\n")
	filter = subprocess.run(['python', os.path.join(directory, 'vhost_parse.py'), \
								'-i', '%s' % (input_fasta), \
								'-o', '%s' % (output_fasta), \
								'-c', '%s' %(vhost_csv), \
								'--taxid'], stdout=subprocess.PIPE, universal_newlines=True).stdout.split("\n")
	try:
		if '' in filter[-1]: del filter[-1]
	except IndexError: # empty list - no taxids are mammalian
		pass
	new_taxids = set(filter)

	return new_taxids

### Output updated db based on output fasta file
def rewrite(output_fasta, new_taxids, db):
	print("(3/3) Outputting updated DB file\n")
	output_csv = str(output_fasta).rsplit('.',1)[0]+".csv"
	#print(output_csv)
	with open(output_csv, 'w+') as out:
		for key,value in db.items():
			#print(value[0])
			if value[0] in new_taxids:
				out.write(str(key)+","+str(value[0])+","+str(value[1])+","+str(value[2]))

### Remove all generated non-results files and directories (including VHost-Classifier)
### TODO: Add flag to keep VHost-Classifier output
def clean():
	print("Cleaning up.\n")
	os.remove(os.path.join(directory, "VHost-Classifier/.vhost_filter.taxid"))
	shutil.rmtree(os.path.join(directory, "VHost-Classifier/vhost_filter"))
	os.remove(os.path.join(directory, "VHost-Classifier/.vhost_filter_mammal.csv"))

### Order of operations
def run():
	parser = create_parser()
	args = parser.parse_args()
	
	db = parse_db(args.input_db_csv)
	vHost(start_dir)
	new_taxids = parse_vHost(args.input_fasta_file, args.output_fasta_file)
	rewrite(args.output_fasta_file, new_taxids, db)
	clean()

if __name__ == '__main__':
### Parameters
	start_dir = os.getcwd()
	script_path = sys.argv[0]
	try:
		directory = script_path.rsplit("/", 1)[0]
	except IndexError:
		directory = os.getcwd() # script in current working directory
	run()
