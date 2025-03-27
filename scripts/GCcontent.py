#!/bin/env python

### Reference: https://gist.github.com/wdecoster/8204dba7e504725e5bb249ca77bb2788
# Takes in a tsv file containing minimally chromosome begin end to calculate the GC content of the targets
# Second argument is the name of the field to take as identifier
# Alternative input is a fasta file
# wdecoster

import sys, os, time
from Bio import Entrez, SeqIO
Entrez.email = "example@email.com"

def getGC(seq):
	return(str(int(round((sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / len(seq)) * 100))))

def getFas(chrID, start, stop):
	record = str(SeqIO.read(Entrez.efetch(db="nucleotide", id=chrID,
		rettype="fasta", strand="1", seq_start=start, seq_stop=stop), "fasta").seq)
	return([getGC(record), str(len(record)), str(record)])

def processFasta(fas):
	print("identifier\tGCcontent\tlength")
	for record in SeqIO.parse(fas, "fasta"):
		print("{}\t{}\t{}".format(record.id, getGC(str(record.seq)), str(len(record.seq))))

def processBed(bedstream, namefield):
	hg19_gidict = {'1': "224589800", '2': "224589811", '3': "224589815", '4': "224589816", '5': "224589817", '6': "224589818",
		'7': "224589819", '8': "224589820", '9': "224589821", '10': "224589801", '11': "224589802", '12': "224589803", '13': "224589804",
		'14': "224589805", '15': "224589806", '16': "224589807", '17': "224589808", '18': "224589809", '19': "224589810", '20': "224589812",
		'21': "224589813", '22': "224589814", 'X': "224589822", 'Y': "224589823", 'M': "17981852"}
	print("identifier\tGCcontent\tlength\tsequence")
	for line in bedstream:
		linelist = line.strip().split()
		print("{}\t{}".format(linelist[namefield], "\t".join(getFas(hg19_gidict[linelist[0].replace('chr','')], linelist[1], linelist[2]))))
		time.sleep(1)

def main():
	if not len(sys.argv) in [2,3]:
		sys.exit("ERROR: Incorrect number of arguments\nEither supply a bed file with header and name of identifier to use, or a fasta file.")
	with open(sys.argv[1], 'r') as input:
		for line in input:
			if line.startswith('#'):
				continue
			else:
				if line.startswith("chromosome\tbegin\tend"):
					try:
						namefield = line.strip().split().index(sys.argv[2]) #Find the index of the field to be used as identifier
						processBed(input, namefield)
						break
					except ValueError:
						sys.exit("Could not find {} in the header".format(sys.argv[2]))
				elif line.startswith(">"):
					processFasta(sys.argv[1])
					break
				else:
					sys.exit("ERROR:\nRequired header not found, expecting minimally 'chromosome' 'begin' 'end' '<identifier>\nAlternatively use a fasta file as input.'")

main()
