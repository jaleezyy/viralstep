#!/bin/env python

### Goal is to use update 'nt' database for use with local BLASTn
### TODO: Download taxdb naming from ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz and remove after script finishes
### Download  'nt' database from ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz

import os, sys, subprocess, shutil
from Bio.Blast.Applications import NcbimakeblastdbCommandline

start_location = os.getcwd()
script_path = os.path.abspath(sys.argv[0])
try:
	update = eval(sys.argv[1]) # update BLASTDB if found
	skip = eval(sys.argv[2])
except NameError: # something that isn't boolean
	update = str(sys.argv[1])
	skip = str(sys.argv[2])
except IndexError: # no argument provided
	update = True
	skip = False

if isinstance(update, bool) is False:
	if "help" in str(update).lower():
		exit("""
	Use this script to update the nt (nucleotide) BLASTDB database.
	
	Run script as:
		python update_blastdb.py <optional_update: True[default]/False> <optional_skip: True/False[default]>
		
	If no arguments provided, then BLASTDB will be updated (in its entirety) by default.
	If update is "False", BLAST+ provided update_blastdb will run to verify complete database.
	If skip is "True" and BLASTDB is found, the script ends.
	If "help" is typed instead, this message will appear.
	
	If you already have a directory called "blastdb" in the directory this script is found, then re-name it if you wish to preserve it.
	""")
	else: # might change to automatically run update = True (accident-prone?)
		exit("""Invalid argument. Re-run with "help" to see options.""")
# else:
	# update = True
	
try:
	directory = script_path.rsplit("/", 1)[0]
except IndexError:
	directory = os.getcwd()

### Check filesystem is properly setup
if os.path.isdir(os.path.join(directory, "blastdb")): 
	print("BLASTDB nt found.")
	os.chdir(os.path.join(directory, "blastdb"))
	if skip: exit("Skipping update of BLASTDB.\n")
else:
	print("Creating new directory for BLASTDB nt.\n")
	os.mkdir(os.path.join(directory, "blastdb"))
	os.chdir(os.path.join(directory, "blastdb"))

### Use NCBI BLASTDB update tool "update_blastdb.pl"
try:
	update_blastdb = subprocess.run(["which", "update_blastdb"], stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.strip("\n")
except subprocess.CalledProcessError:
	print("Error or no update_blastdb.pl found, suggesting BLAST+ (v2.9.0) is not installed.")
	exit()

if update:
	print("\nUpdating BLASTDB.\n")
	subprocess.run(["perl", "%s" %(update_blastdb), "--blastdb_version", "5", "--passive", "--decompress", "--force", "--verbose", "--verbose", "nt", "taxdb"])
else:
	# Note this will technically still "update" the DB - better change below to validate filenames
	print("Verifying BLASTDB.\n")
	subprocess.run(["perl", "%s" %(update_blastdb), "--blastdb_version", "5", "--passive", "--decompress", "--verbose", "--verbose", "nt", "taxdb"])

os.chdir(start_location)
print("Update complete.\n")
