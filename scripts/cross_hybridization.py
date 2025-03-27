#!/bin/env python

from Bio import SeqIO
import argparse
import re
import os, sys, subprocess
import glob
import csv
import datetime, time
import textwrap
from multiprocessing import Pool, Process, Lock, Queue, SimpleQueue, Pipe
from Bio.Blast.Applications import NcbiblastnCommandline
from collections import OrderedDict

### Version 5
### Goal is to filter bait sequences based on BLASTn alignment (with multiprocessing)
### Want to remove BLASTn results linked to Viruses first
### Remaining results are to be used to remove baits with matching elements in the header (i.e. AccID_rng)

### BLASTn should be run with the following: -outfmt '10 std staxids sscinames sblastnames sskingdoms'

### TODO: Add blastn functionality
### TODO: Download taxdb naming from ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz and remove after script finishes
### Can either use remote web BLAST or download  'nt' database from wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz

def create_parser():
	parser = argparse.ArgumentParser(prog="viralstep crosshybaits", description='Filter bait sequences based on BLAST alignment.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file")
	parser.add_argument('-b', '--blast', dest="input_blast_file", default=None, required=False, help="Input BLAST results file; will overwrite automatic BLAST analysis. At minimum the expected file should be .csv with 3 columns: 'Query', 'Reference Hit', 'Sequence similarity %%', 'Alignment length', 'Keep_Term (sskingdom)'")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="Output fasta file")
	parser.add_argument('-t', '--threads', dest="num_threads", required=False, default=1, type=int, help="Number of threads. Default = 1")
	parser.add_argument('-k', '--keep', dest="keep_term", required=False, default="Viruses", help="Keep baits with matching label. Presumed label is sskingdom from BLASTn, but must be found on the last column.")
	parser.add_argument('--validate', dest="validate_baits", action="store_true", required=False, help="Include validation step AFTER off-target filtration to provided bait set (Will take significantly more time to complete).")
	parser.add_argument('--validate-only', dest="only_validate", action="store_true", required=False, help="Only perform validation of given baits from a FASTA file.")
	parser.add_argument('--keep-intermediate', dest="output_pre", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	#parser.add_argument('--updateblastdb', dest="update_blastdb", action="store_true", required=False, help="If no BLAST file provided, and the nt database is found, force an update of the DB (will take additional time). Will supercede '--no-update' flag.")
	#parser.add_argument('--no-update', dest="no_update_blastdb", action="store_true", required=False, help="If no BLAST file provided, and BLASTDB is found, continue with filter.")
	parser.add_argument('--strict', dest="strict_filter", action="store_true", required=False, help="Stricter off-target filter where baits are removed if a single significant off-target is found. Generally faster but you may lose relevant baits.")
	return parser

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
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq
	
def count_line(input):
	num_line = len([1 for line in open(input)])
	return num_line

def generator_file(file):
	with open(file) as r:
		for line in r:
			yield line

def split_file(file, num, threads):
	batch_all = [] # store split input
	if file.endswith(tuple([".fasta", ".fa", "_pre-validate", ".fa_trimmed", ".fasta_trimmed"])):
		record_iter = SeqIO.parse(file, 'fasta')
	else:
		record_iter = iter(generator_file(file))
	if threads <= 0:
		raise ValueError
	for i, batch in enumerate(batch_iterator(record_iter, int(num//threads + (num % threads > 0)))):
		i = i + 1
		if i <= threads:
			batch_all.append(batch)
		else:
			break
	return batch_all
		

### TODO: Output both of these lists to individual output files
### TODO: f.readlines() on keep_list equivalent file to generate same keepage list for validation
### ABOVE REPLACED WITH MULTIPROCESSING.PIPE - still testing
### Set up open files prior to running functions below - done
### Use multiprocessing.async() across a number of processes to speed up execution (refer to syndromic filter) - individual processes instead
### Given all tmp output is here, I don't think I need the Queue to call a writer - order does not matter here - use pipes to retain certain variables
### Use sets to add instead of lists - no duplicates by their very nature
def parse_blast(blast, keep_term, out_remove, queue): # function acts to remove entries containing the sskingdom term (will have to expand)
	for line in blast:
		#print(line)
		if str(keep_term).upper() in str(line).upper().split(",")[-1]:
			print("Keep: "+str(line).split(",")[0])
			queue.put(str(line).split(",")[0]) # send to set object, duplicates not allowed
			#keep_list.append(str(line).split(",")[0]) # used to check bait header
		else:
			print("Check: "+str(line).split(",")[0])
			out_remove.write(line)
			#blast_list.append(line) #full blast line --> to tmp file

def parse_blast_queue(queue, stop_token):
	keepage = set()
	while True:
		blast_keep = queue.get()
		if stop_token in blast_keep:
			catch_blast.put(keepage)
			return
		else:
			keepage.add(blast_keep)

def to_remove(seq, blast): # generates list of sequence headers (assumed AccID_baitrng) that are NOT to be found in the final FASTA file
	removal = set()
	sequences = {}
	for seq_record in SeqIO.parse(seq, 'fasta'):
		sequences[str(seq_record.description)] = str(seq_record.seq).upper()
	
	for line in open(blast):
		query = line.split(",")[0]
		if query in sequences:
			removal.add(query)
		else:
			pass
	return list(removal)

### Convert seq into dict with AccID: Sequence
### Query term is the AccID, so can search dict by AccID
### Read through blast only once
### Compare sequence similarity [2] and alignment length [3] concurrently
### Write directly to file
def off_target(seq, removal_list, blast, out, strict): 
# filters the input FASTA, outputting new FASTA file reduced based on hits found in BLAST output
	count_off_target = 0
	sequences = OrderedDict()
	seq_good = OrderedDict()
	seq_bad = OrderedDict()
	
	for seq_record in SeqIO.parse(seq, 'fasta'):
		sequences[str(seq_record.description)] = str(seq_record.seq).upper()
		seq_good[str(seq_record.description)] = 0
		seq_bad[str(seq_record.description)] = 0
		
	for line in open(blast):
		query = line.split(",")[0]
		print("Query: " + query)
		seq_sim = float(line.split(",")[2])
		align_len = float(line.split(",")[3])
		if query in removal_list:
			if query in line.split(",")[1]: # if hit matches query, ignore
				continue
			if align_len > 50.0 and seq_sim > 80.0:
				seq_bad[query] = seq_bad[query] + 1 # bad off target hit
			else:
				seq_good[query] = seq_good[query] + 1
	
	if strict is True: print("Using strict filter.\n")
	for seq,good,bad in zip(sequences, seq_good, seq_bad):
		assert seq == good
		assert seq == bad
		if strict is True:
			if int(seq_bad[bad]) > 0:
				continue
			else:
				print("Keeping: " + str(seq))
				out.write(">" + str(seq) + "\n" + str(sequences[seq]) + "\n")
				count_off_target = count_off_target + 1
		elif int(seq_good[good]) >= int(seq_bad[bad]): # equal includes 0/0, including virus hits == meaningful hits > unwanted hits
			print("Keeping: " + str(seq))
			out.write(">" + str(seq) + "\n" + str(sequences[seq]) + "\n")
			count_off_target = count_off_target + 1
		else:
			continue
	
	return count_off_target
	
### TODO: Add rejected baits option

### TODO: Heavily requires optimization
### Consider same optimization strategy as off_target function
def validate_check(seq, keep, blast, keep_term, out_file):
	if len(keep) < 1:
		return print("\nNothing to validate.\n")
	else:
		pass
	sequences = OrderedDict()
	seq_score = OrderedDict()
	count_validation = 0
	for seq_record in SeqIO.parse(seq, 'fasta'):
		sequences[str(seq_record.description)] = str(seq_record.seq).upper()
		seq_score[str(seq_record.description)] = 0

	for line in open(blast):
		score = 0 # score metric
		query = line.split(",")[0]
		sskingdom = line.split(",")[-1].upper().strip("\n")
		if query not in sequences: # filtered out from off-targer
			continue
		if query in keep and sskingdom in keep_term.upper():
			print("Scoring %s which has a hit belonging to %s" % (query, line.split(",")[-1].strip("\n")))
			seq_sim = float(line.split(",")[2])
			if seq_sim == float(100.000):
				seq_score[query] = seq_score[query] + 3
			elif float(98) <= seq_sim < float(100):
				seq_score[query] = seq_score[query] + 2
			elif float(95) <= seq_sim < float(98):
				seq_score[query] = seq_score[query] + 1
			elif float(90) <= seq_sim < float(95):
				seq_score[query] = seq_score[query] - 1
			elif float(80) <= seq_sim < float(90):
				seq_score[query] = seq_score[query] - 2
			elif float(70) <= seq_sim < float(80):
				seq_score[query] = seq_score[query] - 3
			elif seq_sim < float(70.000):
				seq_score[query] = seq_score[query] - 4
		else:
			continue
		
		#count_validation = 0
	for acc in seq_score:
		if seq_score[acc] >= 0:
			print("Passed validation: " + str(acc)) 
			out_file.write(">" + str(acc) + "\n" + str(sequences[acc]) + "\n")
			count_validation = count_validation + 1
		
	print("\nValidated baits.\n")
	return count_validation

def blastdb(update, skip, script=None):
### Check if nt database is present and update if asked
	if script is not None:
		directory = script
	else:
		directory = '.'
	print("Checking BLASTDB.\n")
	if update is True:
		subprocess.run(['python', '%s/update_blastdb.py' %(directory), 'True', 'False'])
	elif skip is True:
		subprocess.run(['python', '%s/update_blastdb.py' %(directory), 'False', 'True'])
	else:
		# runs verification
		subprocess.run(['python', '%s/update_blastdb.py' %(directory), 'False', 'False'])

def run_blast(input_fasta, fasta_name, out_dir, threads):
	print("Running BLASTn.\n")
	blastdb = os.path.join(directory, "blastdb")
	if not os.path.isdir(blastdb): 
		exit("Missing BLASTDB! Please run 'viralstep downloadblastdb' or 'update_blastdb.py'!")
	else:
		os.chdir(blastdb)
	blastline = NcbiblastnCommandline(query=input_fasta, db='nt', outfmt="10 std staxids sscinames sblastnames sskingdoms", out='%s/%s.out' % (out_dir, fasta_name), num_threads=threads)
	stdout, stderr = blastline()
	os.chdir(start_dir)
	
def clean(fasta_name, out_dir, out_file, input_blast_file): # remove intermediate files
	if input_blast_file is False:
		if os.path.exists("%s/%s.out" %(out_dir, fasta_name)): os.remove("%s/%s.out" %(out_dir, fasta_name))
	if os.path.exists("%s/%s.out_remove_tmp" %(out_dir, fasta_name)): os.remove("%s/%s.out_remove_tmp" %(out_dir, fasta_name))
	if os.path.exists("%s/%s_pre-validate" %(out_dir, out_file)): os.remove("%s/%s_pre-validate" %(out_dir, out_file))

def run():
### Setup parameters
	processes = []
	parser = create_parser()
	args = parser.parse_args()
	count_off_target = 0
	count_validation = 0
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	try:
		out_dir = os.path.abspath(args.output_fasta_file).rsplit("/", 1)[0]
	except IndexError: # same as cwd
		out_dir = os.getcwd()
	try:
		threads = int(args.num_threads)
	except TypeError:
		print("Invalid thread value. Defaulting to 1.\n")
		threads = 1

### Check if BLAST file provided, otherwise run BLAST with "nt" database
	if args.input_blast_file is None:
		input_blast_file = False
		#blastdb(args.update_blastdb, args.no_update_blastdb)
		run_blast(fasta_file, fasta_name, out_dir, threads)
		blast_file = str("%s/%s.out" %(out_dir, fasta_name)) # replace with os.path
	else:
		blast_file = args.input_blast_file
		input_blast_file = True

	print("Parsing BLAST output.\n")
	
	num_blast = count_line(blast_file)
	batch_blast = split_file(blast_file, num_blast, threads)
	print("Split BLAST input into " + str(len(batch_blast)) + " subunit(s).\n")
	
	num_seq = count_seq(fasta_file)
	
### Parse BLAST output file

	print("Keeping sequences belonging to sskingdom: " + str(args.keep_term) + "\n")
	
	parse_process = Process(target=parse_blast_queue, args=(queue_parse_blast, stop_token))
	parse_process.start()
	
	remove_tmp = open(blast_file+"_remove_tmp", 'w+')
	#keep_tmp = open(blast_file+"_keep_tmp", 'w+')
	for i in range(len(batch_blast)):
		p = Process(target=parse_blast, args=(batch_blast[i],args.keep_term,remove_tmp,queue_parse_blast))
		processes.append(p)
		p.daemon = True
		p.start()
	for p in processes:
		p.join()
	
	queue_parse_blast.put(stop_token)
	keepage = list(catch_blast.get())
	parse_process.join()
	
	blast = blast_file+"_remove_tmp"
	remove_tmp.close()
	processes.clear()

### Check if only validation step is needed - if so, run validate_check only
	
	if args.only_validate is True:
		print("Running validation only. Validating baits.\n")
		out_validation = open(args.output_fasta_file, 'w+')
		count_validation = validate_check(fasta_file, keepage, blast_file, args.keep_term, out_validation)
		count_validation = num_seq - count_validation
		
		if args.output_pre is False:
			print("Cleaning up.\n")
			clean(fasta_name, out_dir, args.output_fasta_file, input_blast_file)
		else:
			print("Intermediate files will NOT be deleted.\n")
			pass
		return count_validation
	else:
		pass

	
### Using BLAST output file (parsed or otherwise), figure out which baits need to be removed	

	print("\nDetermining baits to remove.\n")
	removal = to_remove(fasta_file, blast)
	
### Off-target filter using above list of baits to remove (or ignore when checking each bait)	
	
	print("Filtering baits.\n")
	if args.validate_baits is True:
		out_file = open(args.output_fasta_file+"_pre_validate", 'w+')
		out_validation = open(args.output_fasta_file, 'w+')
	else:
		out_file = open(args.output_fasta_file, 'w+')
	
	count_off_target = off_target(fasta_file, removal, blast, out_file, args.strict_filter)
	count_off_target = num_seq - count_off_target
	
### Cross-validate each bait with the remaining (keep_term) hits (using unparsed BLAST file)	
	
	if args.validate_baits is True:
		print("\nValidating baits.\n")
		### Re-count and recreate sequence split (with filtered sequences)
		
		fasta_file = args.output_fasta_file+"_pre_validate"
		num_seq_2 = count_seq(fasta_file)
		
		count_validation = validate_check(fasta_file, keepage, blast_file, args.keep_term, out_validation)
		count_validation = num_seq_2 - count_validation
	else:
		print("\nSkipping validation step.\n")

### Clean up intermediate files unless otherwise asked not to
### TODO: Move tmp files to central directory (at cwd)
	if args.output_pre is False:
		print("\nCleaning up.\n")
		clean(fasta_name, out_dir, args.output_fasta_file, input_blast_file)
	else:
		print("Intermediate files will NOT be deleted.\n")
		pass

	return (count_off_target + count_validation)

if __name__ == '__main__':
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	script_path = os.path.abspath(sys.argv[0])
	queue_parse_blast = SimpleQueue() # update [keepage]
	catch_blast = SimpleQueue() # catch final list (may replace with single queue)
	stop_token = "STOP"
	start_dir = os.getcwd()
	try:
		directory = script_path.rsplit("/", 1)[0]
	except IndexError:
		directory = os.getcwd() # script in current working directory
	num_baits = run()
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	print("\nRemoved " + str(num_baits) + " bait(s).\n")
	print("Start: " + start)
	print("Finish: " + finish)
	exit()
