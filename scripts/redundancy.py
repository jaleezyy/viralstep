#!/bin/env python

### Goal is to create my own version of the greedy algorithm based on BLAST local alignment
### Run BLASTn on bait set itself
### Isolate a set of bait headers (SeqRecord) - check memory usage
### Parse BLAST file relative to bait header set
### Will integrate BLASTDB generation, alignment, and analysis in one script (future)

import sys, os, subprocess, shutil
import argparse
from Bio import SeqIO
import datetime, time
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

def create_parser():
	parser = argparse.ArgumentParser(prog="viralstep redundbaits", description='Removes redundant baits based on self-BLAST alignment.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file")
	parser.add_argument('-l', '--length', dest="align_length", required=False, help="Alignment length of hit based on BLASTn. '-s' (sequence similarity) required. Incompatible with '-m' (sequence identity).")
	parser.add_argument('-s', '--similarity', dest="align_sim", required=False, help="Sequence similarity of hit based on BLASTn. '-l' (alignment length) required. Incompatible with '-m' (sequence identity).")
	parser.add_argument('-m', '--sequence-identity', dest="sequence_ident", required=False, help="Percent sequence identity. Script will estimate sequence identity between baits and run filter based on given input percent matching nucleotides. Not compatible with '-s' (sequence similarity) and/or '-l' (sequence length).")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="Output fasta file")
	parser.add_argument('-t', '--threads', dest="num_threads", required=False, default=1, type=int, help="Number of threads for BLASTn. Default = 1")
	parser.add_argument('--keep-intermediate', dest="output_pre", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	return parser

def count_seq(input):
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq
	
def makeselfblastdb(script_dir, fasta_file, out_dir):
	fasta_file = os.path.abspath(fasta_file)
	db_path = os.path.join(script_dir, "self_blastdb")
	if os.path.isdir(db_path):
		print("Updating existing self-BLAST DB\n")
		shutil.rmtree(db_path) # delete previous self-blast
		os.mkdir(db_path)
		os.chdir(db_path)
	else:
		os.mkdir(db_path)
		os.chdir(db_path)
	print("Generating BLASTDB of provided input FASTA file.\n")
	dbline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=fasta_file, out="self-blast")
	stdout, stderr = dbline()
	
def selfblastn(fasta_file, fasta_name, script_dir, out_dir, threads):
	### Running from scripts_dir/self_blastdb
	# if os.getcwd() != os.path.abspath(os.path.join(script_dir, "self_blastdb")):
		# os.chdir(os.path.abspath(os.path.join(script_dir, "self_blastdb")))
	#fasta_path = os.path.abspath(fasta_file)
	blast = os.path.join(out_dir, fasta_name+".out")
	print("Running BLASTn\n")
	blastline = NcbiblastnCommandline(query=fasta_file, db='self-blast', outfmt='10 std staxids sscinames sblastnames sskingdoms', out=blast, num_threads=threads)
	stdout, stderr = blastline()
	os.chdir(start_dir)
	return blast

def greedy(bait_set, bait_hits, align_sim, align_length): # uses seq_sim and align_length
	for bait in bait_hits:
		if bait not in bait_set: # if already removed, don't analyse
			continue
		else:
			hits_list = bait_hits[bait]
		for i in hits_list:
			length = float(str(i).split(",")[2])
			similarity = float(str(i).split(",")[1])
			if similarity >= float(align_sim) and length >= float(align_length):
				try:
					del bait_set[i.split(",")[0]]
				except KeyError:
					print("Bait " + str(i.split(",")[0]) + " not found, skipping.")
					continue
			else:
				continue
	return bait_set
	
def greedy_ident(bait_set, bait_hits, sequence_iden, bait_length): # calculates percent identical nucleotides
	for bait in bait_hits:
		if bait not in bait_set: # if already removed, don't analyse
			continue
		else:
			hits_list = bait_hits[bait]
		for i in hits_list:
			length = float(str(i).split(",")[2])
			similarity = float(str(i).split(",")[1])/100
			percent_identical = ((similarity*length)/bait_length)*100
			if percent_identical >= float(sequence_iden):
				try:
					del bait_set[i.split(",")[0]]
				except KeyError:
					print("Bait " + str(i.split(",")[0]) + " not found, skipping.")
					continue
			else:
				continue
	return bait_set

def clean(blast): # remove intermediate files
	if os.path.exists(blast): os.remove(blast)
	
def run():
	### Parser
	parser = create_parser()
	args = parser.parse_args()

	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	num_seq = count_seq(fasta_file)
	ident_calc = None # check whether filter calculates sequence identity
	try:
		out_dir = os.path.abspath(args.output_fasta_file).rsplit("/", 1)[0]
	except IndexError:
		out_dir = os.getcwd()
	threads = int(args.num_threads)
	if (args.align_sim is not None or args.align_length is not None):
		if args.sequence_ident is not None:
			exit("'-m' (sequence_ident) is not compatible with '-l' (align_length) and '-s' (align_sim). Refer to '-h/--help' for more information.")
		try:
			align_sim = float(args.align_sim)
			align_length = float(args.align_length)
			ident_calc = False
		except (IndexError, ValueError, TypeError):
			exit("'-l' (align_length) and '-s' (align_sim) are both required Refer to '-h/--help' for more information.")
	elif args.sequence_ident is not None:
		if (args.align_sim is not None or args.align_length is not None):
			exit("'-m' (sequence_ident) is not compatabile with '-l' (align_length) and '-s' (align_sim). Refer to '-h/--help' for more information.")
		try:
			sequence_ident = float(args.sequence_ident)
			ident_calc = True
		except (ValueError, TypeError):
			exit("Invalid value for '-m' (sequence_ident). Refer to '-h/--help' for more information")
	else:
		exit("Missing parameters. Refer to '-h/--help' for more information.")
		
	### Run BLASTn against input fasta file (self-BLAST)
	makeselfblastdb(script_path, fasta_file, fasta_name)
	blast = selfblastn(fasta_file, fasta_name, script_path, out_dir, threads)
	
	bait_set = {}
	bait_length = 0
	for seq_record in SeqIO.parse(fasta_file, 'fasta'): # extract all baits - subject to removal)
		bait_set[seq_record.description] = seq_record.seq.upper()
		if (bait_length == 0 or bait_length > len(str(seq_record.seq))):
			bait_length = len(str(seq_record.seq)) # final bait length should be lowest common denominator

	print("Processed FASTA file.\n")	
	#print("Bait Set Before: \n")
	#print(bait_set)
	#print("\n")

	print("Running optimized BinB4Greedy Algorithm!\n")
	bait_hits = defaultdict(list) # Query: [hits]
	#analysed = set()
	for line in open(blast): 
		if str(line).split(",")[0].split("_")[0] in str(line).split(",")[1]:
			continue
		else:
			print("Bait " + str(line).split(",")[0] + " hits Bait " + str(line).split(",")[1] + "!")
			bait_hits[str(line).split(",")[0]].append(str(line).split(",",1)[1])
	
	if ident_calc:
		bait_set = greedy_ident(bait_set, bait_hits, sequence_ident, bait_length)
	else:
		bait_set = greedy(bait_set, bait_hits, align_sim, align_length)

	count_baits = 0
	with open(args.output_fasta_file, 'w+') as out:
		for bait in bait_set:
			out.write(">" + str(bait) + "\n" + str(bait_set[bait]) + "\n")
			count_baits = count_baits + 1
	count_baits = num_seq - count_baits
	
	### Clean up intermediate files unless otherwise asked not to
	if args.output_pre is False:
		print("\nCleaning up.\n")
		clean(blast)
	else:
		print("\nIntermediate files will NOT be deleted.\n")
		pass
	
	print("Finished removing redundant probes.\n")
	
	return count_baits
			
				


if __name__ == '__main__':
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	try:
		script_path = sys.argv[0].rsplit("/",1)[0]
	except IndexError:
		script_path = os.getcwd() # script in current working directory
	start_dir = os.getcwd()
	num_baits = run()
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	print("\nRemoved " + str(num_baits) + " bait(s).\n")
	print("Start: " + start)
	print("Finish: " + finish)
	exit()
		
		
	
	
	
