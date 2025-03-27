#!/bin/env python

import os 
import subprocess
import sys
from Bio import SeqIO
#from Bio.SeqUtils.CheckSum import seguid

### REQUIRES OLIGOARRAYAUX 3.8
# melt_parse.sh executes melt.pl with subsequence parsing (MUST BE IN SAME PATH AS THIS PYTHON SCRIPT)
# After running melt.pl determine which probes are bad relative to melting temperature (Tm)
# Takes 2 files: probes and melting point (mp) information
# Will compare each value in mp and if greater than specified threshold then it will take the header + sequence from probes to new filtered file

HELP="""
Melting temperature filter for a given bait set.
Requires OligoArrayAux 3.8 and ViralSTEP melt_parse_baits.sh.

Arguments are as followed:
	melting_temperature.py <bait_set_fasta> <min_temp (default = 65)> <max_temp (default = 95)>
	
Note:
    - If list of melting temperature with filename: Tm_<fasta_filename.fasta>.txt found, then melt_parse_baits.sh will be skipped for time
    - Temperatures are in deg C.
    - Output will be found in same directory as input FASTA file.

"""

try:
	if 'help' in sys.argv: exit(HELP)
	script_path = os.path.join(os.path.abspath(sys.argv[0]).rsplit("/",1)[0], "melt_parse_baits.sh") # script directory
	userParameters = sys.argv[1] # input FASTA file
	low_temp = sys.argv[2] # input minimum temperature
	high_temp = sys.argv[3] # input maximum temperature
except IndexError:
	try:
		if 'help' in sys.argv: exit(HELP)
		script_path = os.path.join(os.path.abspath(sys.argv[0]).rsplit("/",1)[0], "melt_parse_baits.sh") # script directory
		userParameters = sys.argv[1] # input FASTA file
		low_temp = 65 # input minimum temperature (hybridization)
		high_temp = 95 # input maximum temperature (denaturation)
		#print("Missing parameters! Attempting to run with default temperature range: 65-95 deg C.\n")
	except IndexError:
		exit(HELP)

if not os.path.exists(userParameters): exit("Cannot find input FASTA file!")

try:
	fastafile = userParameters.rsplit("/",1)[1]
	fastapath = userParameters.rsplit("/",1)[0]
	#fastafile = os.path.basename(userParameters)
	#fastapath = os.path.dirname(userParameters)
except IndexError: # current working directory
	fastafile = userParameters
	fastapath = os.getcwd()

# if statement determining whether Tm data already exists as proper name --> otherwise run melt_parse.sh
# allow for user to skip melt_parse.sh if already done
if os.path.exists("%s/Tm_%s.txt" %(fastapath, fastafile)):
	print("Melting temperature data exists for original FASTA file, skipping melt_parsh.sh" + "\n")
	Tm_userParameters = os.path.join(fastapath,"Tm_"+fastafile+".txt")
	clean_required = False
else:
	print("Running melt.pl on given probe set")
	# subprocess.run("~/jalees/Scripts/melt_parse.sh -s %s" %(userParameters), check=True, shell=True)
	subprocess.run(["%s" %(script_path), "-s", "%s" %(userParameters), "-o", "%s/Tm_%s.txt" %(fastapath, fastafile)])
	print("Generated list of melting temperatures" + "\n") # Tm_userParameters.txt
	Tm_userParameters = os.path.join(fastapath,"Tm_"+fastafile+".txt") 
	clean_required = True


print("Original dataset:")
print(userParameters + "\n") # original fasta file with probes


with open(Tm_userParameters, "r") as f:
	lines = f.read().splitlines()
	

print("Melting temperature data:")
print(Tm_userParameters + "\n")
#print(len(lines)) # check length of the list from Tm data which should match length of original fasta file
#print(lines[0]) # check if first value shows

def melt_filter_range(fastapath, fastafile, low_temp, high_temp):
	# Create our hash table to add the sequences
	sequences={} # set object to eliminate duplicates
	no_seq={} 
	i = 1 # first element in list to examine

	if isinstance(low_temp, str) == True:
		low_temp = float(low_temp)
		print("Minimum temperature: " + str(low_temp))
	else:
		print("Error! Invalid minimum temperature!")
		exit()
		
	#print("Maximum temperture:")
	#high_temp = input()
	
	if isinstance(high_temp, str) == True:
		high_temp = float(high_temp)
		print("Maximum temperature: " + str(high_temp))
	else:
		print("Error! Invalid maximum temperature!")
		exit()
	
	header = []
	sequence = []
	no_header = []
	no_sequence = []

	# Using the Biopython fasta parse we can read our fasta input
	for seq_record in SeqIO.parse(os.path.join(fastapath,fastafile), "fasta"):
		#sequence = str(seq_record.seq).upper()
		if float(lines[i]) >= low_temp and float(lines[i]) <= high_temp:
			header.append(seq_record.description)
			sequence.append(seq_record.seq)
			#sequences[sequence] = seq_record.description
			i = i + 1
		else:
			no_header.append(seq_record.description)
			no_sequence.append(seq_record.seq)
			#no_seq[sequence] = seq_record.description
			#next
			i = i + 1

	# Create output file 
	output_pass = os.path.join(fastapath,"TmFiltered_"+fastafile)
	output_fail = os.path.join(fastapath,"rejected_TmFiltered_"+fastafile)
	
	# Create a file in the same directory where script was run
	with open(output_pass, "w+") as output_file:
		# Just read the hash table and write on the file as a fasta format
		for seq,head in zip(sequence,header):
			output_file.write(">" + str(head) + "\n" + str(seq) + "\n")
	with open(output_fail, "w+") as output_file:
		for seq,head in zip(no_sequence,no_header):
			output_file.write(">" + str(head) + "\n" + str(seq) + "\n")
	
	print("\nRemoved bad probes!\nOutput file: %s\nUnwanted sequences found in output file: %s" %(output_pass, output_fail))
	
def clean(fastapath, fastafile, start_dir, Tm_userParameters):
	### Remove intermediate files from melt.pl
	file = os.path.join(fastapath, fastafile)
	try: # if linking to older melt.pl bash script
		os.remove("%s.65.ext" %(file))
		os.remove("%s.65.plot" %(file))
		os.remove("%s.ct" %(file))
		os.remove("%s.dG" %(file))
		os.remove("%s.run" %(file))
	except (OSError, FileNotFoundError): # if using newer symlink melt.pl bash script
		os.remove(os.path.join(start_dir, "%s.symlink.65.ext" %(fastafile)))
		os.remove(os.path.join(start_dir, "%s.symlink.65.plot" %(fastafile)))
		os.remove(os.path.join(start_dir, "%s.symlink.ct" %(fastafile)))
		os.remove(os.path.join(start_dir, "%s.symlink.dG" %(fastafile)))
		os.remove(os.path.join(start_dir, "%s.symlink.run" %(fastafile)))
	print("\nUsing %s, you can re-run the filter for a different temperature range.\nEnsure both the original FASTA and melting temperatures are in the same directory if running this script directly, otherwise use 'viralstep tmbaits --tmlist %s." %(Tm_userParameters, Tm_userParameters))


start_dir = os.getcwd()
melt_filter_range(fastapath, fastafile, low_temp, high_temp)
if clean_required: clean(fastapath, fastafile, start_dir, Tm_userParameters)

