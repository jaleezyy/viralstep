#!/bin/env python

### Goal is to recreate classification pipeline through the creation of a central file
### Central file contains all known virus sequence accessionID, taxIDs, virus families

### Use two dictionary objects: one with taxID to tuple of accID; one with taxID to virus family
### First generate a list of taxID, virus families for the input FASTA file - (tuples?)
### Compare index for virus family in pre-defined groups/families
### Pull corresponding seq_record and output directly to file (rather than to another list)

### This version aims to incorporate multiprocessing to improve efficiency

import argparse
from Bio import SeqIO
import os, sys
import datetime
from collections import defaultdict
from multiprocessing import Process, Queue

def create_parser():
	parser = argparse.ArgumentParser(prog="viralstep classify", formatter_class=argparse.RawDescriptionHelpFormatter, description='Obtain virus families for a given FASTA file in format.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input fasta file")
	parser.add_argument('-l', '--list', dest="input_family_file", required=True, help="Input reference list of virus families identified for each reference sequence. Columns should be: AccessionID, TaxID, Baltimore classificaton, Virus family.")
	parser.add_argument('-o', '--output', dest="output_file", required=True, help="Output CSV filename. WIll have 'family_' attached as prefix.")
	#parser.add_argument ('-f', '--filter', dest="filter_name", required=True, help="Specify virus family or group of virus families (See broad groups below)")
	#parser.add_argument('-v', '--rejected', dest="not_filtered", action="store_true", required=False, help="Obtain remaining sequences (wash-through) rejected by filter.")
	parser.add_argument('-t', '--threads', dest="threads", required=False, help="Number of threads. Will attempt to assign individual threads for each output file. Default = 1.")
	parser.add_argument('--taxid', dest="receive_taxid", action="store_true", required=False, help="Obtain taxID list as additional output")
	#parser.add_argument('--lineage', dest="receive_lineage", action="store_true", required=False, help="Obtain filtered virus family list, derived from lineage data, as additional output")
	parser.add_argument('--log', dest="output_log", action="store_true", required=False, help="Obtain log file output")
	return parser

def load_reference(args):
	accID2taxID = defaultdict(list)
	taxID2family = {} # normal dictionary object since one taxID = one family
	family2baltimore = {} # family: baltimore classification (i.e. dsDNA, +ssRNA, etc.)
	family_list = args.input_family_file
	
	with open(family_list) as f:
		for line in f:
			accID2taxID[line.split(",")[1]].append(line.split(",")[0])
			taxID2family[line.split(",")[1]] = line.split(",")[-1]
			family2baltimore[line.strip("\n").split(",")[-1]] = line.split(",")[-2]
	
	return accID2taxID, taxID2family, family2baltimore

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
	num = len([1 for line in open(input) if line.startswith(">")])
	return num

def split_fasta(fasta_file, num, threads):
	batch_all = [] # store split input
	record_iter = SeqIO.parse(fasta_file,"fasta")
	if threads <= 0:
		raise ValueError("Invalid number of threads")
	for i, batch in enumerate(batch_iterator(record_iter, int(num//threads + (num % threads > 0)))):
		i = i + 1
		if i <= threads:
			batch_all.append(batch)
		else:
			break
	return batch_all
			

def writer_taxID(output, queue_taxID, stop, args):
	# if "." not in output: 
		# output = output + ".txt"
	# elif ".fasta" in output:
		# output = output.rsplit(".",1)[1] + ".txt" 
	# else:
		# pass
	with open(output, "w+") as out_tax:
		while True:	
			line_taxID = queue_taxID.get()
			if stop in line_taxID:
				return
			else:
				out_tax.write(line_taxID)

def writer_lineage(output, queue_lineage, stop, args): # generate output through queues
	# if "." not in output: 
		# output = output + ".txt"
	# elif ".fasta" in output:
		# output = output.rsplit(".",1)[0] + ".txt" # assuming name.fasta = multiple '.' will skew this
	# else:
		# pass
### Output core file (lineage list)		
	with open(output, "w+") as out_lineage:
		while True:
			line_lineage = queue_lineage.get()
			if stop in line_lineage:
				return
			else:
				out_lineage.write(line_lineage)

def writer_alt(output_lineage, output_taxID, queue_lineage, queue_taxID, stop, args): # generate output through queues
	# output_file=output
	# if "." not in output_file: # ensure final output is a .fasta file
		# output_file = output_file + ".txt"
	# elif ".fasta" in output_file:
		# output_file = output_file.rsplit(".",1)[0] + ".txt" 
	# else:
		# pass
### Output core file (lineage list)		
	with open(output_lineage, "w+") as out_lineage:
		if args.receive_taxid is True:
			writer_taxID(output_taxID, queue_taxID, stop, args)
		else:
			pass
		while True:
			line_lineage = queue_lineage.get()
			if stop in line_lineage:
				return
			else:
				out_lineage.write(line_lineage)
		
def classify(batch, args, threads, out_lineage, out_taxID):
### Parameter setup
	accID2taxID, taxID2family, family2baltimore = load_reference(args)
	input_file = batch
	output_file=args.output_file
	if ".csv" not in output_file: # ensure final output is a .csv file
		output_file = output_file + ".csv"
	elif ".fasta" in output_file:
		output_file = output_file.split(".")[0] + ".csv"
	else:
		pass
		
### Obtain and filter sequences
### Check accessionID and search in accID2taxID for corresponding taxID
### Use obtained taxID to search key in taxID2family
### assign family as search_term
### compare and if it matches, take seq_record
### Output directly to file
	

### Single thread file setup 

	if threads < 2:
		out_lineage = open(out_lineage, "w+")
		if args.receive_taxid is True:
			out_tax = open(out_taxID, "w+")
		else:
			pass	
	else:
		pass
	
	search_term = None
	
	for seq_record in input_file:
		for x,y in accID2taxID.items():
			# isolates AccessionID from any sequence header (including post BaitsTools)
			# If interpreting bait header
			if seq_record.description.split(" ")[0].rsplit("_",1)[0] in y: 
				taxID = x
				search_term = str(taxID2family[taxID]).split("\n")[0] # virus family
			# If interpreting regular accession ID
			# Opposite true if db file contains bait headers
			elif seq_record.description.split(" ")[0] in y:
				taxID = x
				search_term = str(taxID2family[taxID]).split("\n")[0] # virus family
			else:
				continue

		if search_term is None:
			print("Unknown")
			if threads < 2:
				out_lineage.write(str(seq_record.id) + ",0,Unknown,Unknown," + str(seq_record.description.split(" ")[0].rsplit("_",1)[0]) + "\n")
			else:
				queue_lineage.put(str(seq_record.id) + ",0,Unknown,Unknown," + str(seq_record.description.split(" ")[0].rsplit("_",1)[0]) + "\n")
### Optional output
			if args.receive_taxid is True:
				if threads < 2:
					out_tax.write("0" + "\n")
				else:
					queue_taxID.put("0" + "\n")
		else:
			print(search_term)
			baltimore_term = family2baltimore[search_term]

			if threads < 2:
				#out_lineage.write(x + "," + search_term + "\n")
				out_lineage.write(str(seq_record.id) + "," + taxID + "," + baltimore_term + "," + search_term + "\n")
			else:
				#queue_lineage.put(x + "," + search_term + "\n")
				queue_lineage.put(str(seq_record.id) + "," + taxID + "," + baltimore_term + "," + search_term + "\n")
### Optional output
			if args.receive_taxid is True:
				if threads < 2:
					out_tax.write(taxID + "\n")
				else:
					queue_taxID.put(taxID + "\n")

### Close optional output files
	if threads < 2:
		out_lineage.close()
		if args.receive_taxid is True:
			out_tax.close()
		else:
			pass
	else:
		pass
				
	
### Check output file
### QC on output FASTA file i.e. check that values are present - no empty file(s)
### Output before and after numbers of sequences

def log_file(args, start, finish):
	with open(args.output_file + ".log", "w+") as output_log:
		output_log.write("Start: " + start + "\n" + "Finish: " + finish + "\n")
	print("\n" + "Start: " + start)
	print("Finish: " + finish)
	print("Filtered sequences based on virus families. Output file: " + "{family,taxID}_" + args.output_file)

	
def run():
	processes = []
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	parser = create_parser()
	args = parser.parse_args()
	input = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(input))
	try:
		out_dir = os.path.abspath(args.output_file).rsplit("/", 1)[0]
	except IndexError: # same as cwd
		out_dir = os.getcwd()
	output_name = str(os.path.basename(args.output_file))
	if output_name.startswith(("family_", "taxID_")):
		output_name = output_name.split("_",1)[1]
	lineage_path = os.path.join(out_dir,"family_%s" %(output_name))
	taxID_path = None
	
	num_seq = count_seq(input)
	expected_output = 1
	if args.receive_taxid is True:
		expected_output = expected_output + 1
		taxID_path = os.path.join(out_dir, "taxID_%s" %(output_name))
	if args.threads is not None:
		try: # multi-writer process	
			if int(args.threads) > int(num_seq): # if true then num_seq << threads (we don't need all threads as a result)
				if int(num_seq-expected_output) <= 0: # if true then num_seq <= expected_output
					raise Exception("Number of sequences << number of threads assigned")
				else:
					args.threads = int(num_seq) # set number of threads equal to num_seq 
					batch = split_fasta(input, num_seq, int(args.threads)-expected_output)
					print("Split input FASTA across %s thread(s), based on number of sequences, dedicating %s thread(s) for writing" % (int(args.threads)-expected_output, expected_output))
			else:
				batch = split_fasta(input, num_seq, int(args.threads)-expected_output) # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
				print("Split input FASTA across %s thread(s), dedicating %s thread(s) for writing" % (int(args.threads)-expected_output, expected_output))
			
###	Writer protocol		
			writer_process = Process(target=writer_lineage, args = (lineage_path, queue_lineage, stop_token, args))
			writer_process.start()
			if args.receive_taxid is True:
				taxID_write = Process(target=writer_taxID, args = (taxID_path, queue_taxID, stop_token, args))
				taxID_write.start()
			else:
				pass
### Process class			
			for i in range(len(batch)): # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
				p = Process(target=classify, args=(batch[i],args, int(args.threads), lineage_path, taxID_path))
				processes.append(p)
				p.deamon = True
				p.start()
			for p in processes:
				p.join()			
			
			queue_taxID.put(stop_token)
			queue_lineage.put(stop_token)
			
			writer_process.join()
			
			if args.receive_taxid is True:
				taxID_write.join()
			else:
				pass
		
		except ValueError: # just give up
			print("Not enough processes for simultaneous writing of output files. Running single process!")
			batch = split_fasta(input, num_seq, int(1))
			classify(batch[0], args, int(1), lineage_path, taxID_path)
		
		except Exception: # one writer process
			try:
				if int(num_seq-1) <= 0:
					raise ValueError("Only one sequence (and therefore) thread available.")
				else:
					args.threads = int(num_seq) 
					batch = split_fasta(input, num_seq, int(args.threads)-1)
					#batch = split_fasta(input, num_seq, int(args.threads)-1) # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
					print("Split input FASTA across %s thread(s), based on number of sequences, with 1 thread dedicated to writing." % (int(args.threads)-1))
				
				writer_process = Process(target=writer_alt, args = (lineage_path, taxID_path, queue_lineage, queue_taxID, stop_token, args))
				writer_process.start()
			
				for i in range(len(batch)): # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
					p = Process(target=classify, args=(batch[i],args, int(args.threads), lineage_path, taxID_path))
					processes.append(p)
					p.deamon = True
					p.start()
				for p in processes:
					p.join()			
			
				queue_taxID.put(stop_token)
				queue_lineage.put(stop_token)
			
				writer_process.join()
				
			except ValueError: # just give up
				print("Cannot split. Running single process!")
				batch = split_fasta(input, num_seq, int(1))
				classify(batch[0], args, int(1), lineage_path, taxID_path)
	
	else: # when threads flag not specified
		print("No thread number specified. Running single process!")
		batch = split_fasta(input, num_seq, int(1))
		classify(batch[0], args, int(1), lineage_path, taxID_path)
			
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	
	if args.output_log is True:
		log_file(args, start, finish)
	else:
		print("\n" + "Start: " + start)
		print("Finish: " + finish)
		print("Updated classification. Output file: " + "{family,taxID}_%s" %(output_name))

if __name__ == '__main__':
	queue_taxID = Queue()
	queue_lineage = Queue()
	stop_token = "Stop!"
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	try:
		directory = script_path.rsplit("/", 1)[0]
	except IndexError:
		directory = os.getcwd() # script in current working directory
	run()
	exit("\n" + "Exited at: " + str(datetime.datetime.utcnow()).split('.')[0])