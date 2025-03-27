#!/bin/env python

### Aim is to take V-Host Classifier output and filter corresponding sequences -- version 4
### Difference between version 3 and 4 is the lack of stdout except list of taxids (from --taxid)

from Bio import SeqIO
#import re
#import sys
import os, argparse
import csv
#import datetime

def create_parser():
	parser = argparse.ArgumentParser(prog='viralstep vhost_parse', description='Filter sequences based on infection host. Run after VHost-Classifier.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="input fasta file")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="output fasta file")
	parser.add_argument('-c', '--csv', dest="input_vhost_csv", required=True, help="input .csv file from VHost-Classifier")
	parser.add_argument('--taxid', dest="receive_taxID", action="store_true", required=False, help="receive list of taxIDs from filtered list")
	return parser
	
def v_host_classify(args):
	#start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	fasta_file = args.input_fasta_file
	if args.output_fasta_file.endswith((".fasta", ".fa")):
		out_file = args.output_fasta_file
	else:
		out_file = args.output_fasta_file.rsplit('.')[0]+".fasta"
	vhost_file= args.input_vhost_csv
	
	seq = []
	seq_narrow=[]
	vhost={}
	vhost_position=[]

	
	for seq_record in SeqIO.parse(fasta_file, "fasta"): # create list of sequences and corresponding header
		seq.append(seq_record)
	
	with open(vhost_file, newline='') as csvfile:
		reader=csv.reader(csvfile, delimiter=",")
		for row in reader:
			vhost[row[1]] = row[0] # seq position: taxID
			
			#vhost_ID.append(row[0])
			vhost_position.append(row[1]) # list of seq positions
			
	
		total = len(vhost_position)
		count = 0
	#	print(str(count) + "/" + str(total))
		vhost_position.sort() # sort prior to indexing

### Create list of filtered sequences
		
	for i in vhost_position:
		indexing=int(i)
		
		seq_narrow.append(seq[indexing])
		count = count + 1
	#	print(str(count) + "/" + str(total))

### Create output file
	# Create a file in the same directory where script was run
	with open(out_file, "w+") as output_file:
		# Just read the hash table and write on the file as a fasta format
		for i in range(0,len(seq_narrow)): #use indexing position 
			indexing=int(i)
			output_file.write(">" + str(seq_narrow[indexing].description) + "\n" + str(seq_narrow[indexing].seq) + "\n")
			
	#finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	
	# with open(out_file + ".log", "w+") as output_log:
		# output_log.write("Start: " + start + "\n" + "Finish: " + finish + "\n")
		
	#print("Start: " + start)
	#print("Finish: " + finish)
	#print("Filtered sequences based on VHost-Classifier. Output file: " + out_file)

### return list of taxIDs
	if args.receive_taxID is True:
	#	print("Extracting filtered taxIDs.")
		sort_position=[] # only when flag is triggered
		for i in vhost:
			sort_position.append(i)
		sort_position.sort()
		for tax in sort_position:
			print(str(vhost[tax]))
		# with open("taxID_" + out_file.split(".fasta")[0] + ".txt", "w+") as output_taxid:
			# for tax in sort_position:
				# output_taxid.write(str(vhost[tax]) + "\n")
	else:
		pass

	
def run():
	parser = create_parser()
	args = parser.parse_args()
	v_host_classify(args)

if __name__ == '__main__':
	run()