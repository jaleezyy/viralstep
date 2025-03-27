#!/bin/env python

### Goal is to recreate classification pipeline through the creation of a central file
### Central file contains all known virus sequence accessionID, taxIDs, virus families
### Filter by family or syndrome
### Use two dictionary objects: one with taxID to tuple of accID; one with taxID to virus family
### First generate a list of taxID, virus families for the input FASTA file - (tuples?)
### Compare index for virus family in pre-defined groups/families
### Pull corresponding seq_record and output directly to file (rather than to another list)

### This version aims to incorporate multiprocessing to improve efficiency

### TODO: Add more flexible -f where a list of virus families (comma-separated) can be submitted


import argparse
from Bio import SeqIO
import os, sys
import datetime
from collections import defaultdict
from multiprocessing import Process, Queue

def create_parser():
	USAGE='''%(prog)s [<args>]

	Broad groups (refer to -f or --filter filter_name argument)
	-----------------------------------------------------------------------------------
		
	'human' (filtering for human-relevant viruses) includes the following virus families:
		- Adenoviridae    - Hantaviridae      - Peribunyaviridae    - Smacoviridae
		- Anelloviridae   - Hepadnaviridae    - Phenuiviridae       - Tobaniviridae
		- Arenaviridae    - Hepeviridae	      - Picobirnaviridae    - Togaviridae
		- Astroviridae    - Herpesviridae     - Picornaviridae
		- Bornaviridae    - Nairoviridae      - Pneumoviridae	
		- Caliciviridae   - Orthomyxoviridae  - Polyomaviridae
		- Circoviridae    - Papillomaviridae  - Poxviridae
		- Coronaviridae   - Paramyxoviridae   - Reoviridae
		- Filoviridae     - Partitiviridae    - Retroviridae
		- Flaviviridae    - Parvoviridae      - Rhabdoviridae
		
	'gastro' includes:
		- Adenoviridae    - Hepeviridae
		- Anelloviridae   - Nairoviridae
		- Arenaviridae    - Orthomyxoviridae
		- Astroviridae    - Partitiviridae
		- Caliciviridae   - Parvoviridae
		- Circoviridae    - Picobirnaviridae
		- Coronaviridae   - Picornaviridae
		- Filoviridae     - Polyomaviridae
		- Hantaviridae    - Reoviridae
		- Hepadnaviridae  - Tobaniviridae
		
	'blood' includes:
		- Adenoviridae    - Nairoviridae      - Reoviridae
		- Anelloviridae   - Papillomaviridae  - Retroviridae
		- Arenaviridae    - Paramyxoviridae   - Rhabdoviridae
		- Bornaviridae    - Parvoviridae      - Togaviridae
		- Filoviridae     - Peribunyaviridae
		- Flaviviridae    - Phenuiviridae
		- Hantaviridae    - Picornaviridae
		- Hepadnaviridae  - Pneumoviridae
		- Hepeviridae     - Polyomaviridae
		- Herpesviridae   - Poxviridae
			
	'respir' includes:
		- Adenoviridae    - Orthomyxoviridae
		- Anelloviridae   - Paramyxoviridae
		- Arenaviridae    - Parvoviridae
		- Bornaviridae    - Peribunyaviridae
		- Coronaviridae   - Phenuiviridae	
		- Filoviridae     - Picornaviridae
		- Flaviviridae    - Pneumoviridae
		- Hantaviridae    - Polyomaviridae
		- Herpesviridae   - Reoviridae
		- Nairoviridae    - Togaviridae
			
	'csf' (Cerebrial Spinal Fluid) includes:
		- Adenoviridae    - Nairoviridae      - Togaviridae
		- Anelloviridae   - Orthomyxoviridae
		- Arenaviridae    - Paramyxoviridae
		- Bornaviridae    - Parvoviridae
		- Coronaviridae   - Peribunyaviridae
		- Filoviridae     - Phenuiviridae
		- Flaviviridae    - Picornaviridae
		- Hantaviridae    - Pneumoviridae
		- Hepeviridae     - Polyomaviridae
		- Herpesviridae   - Rhabdoviridae
			
	'urine' includes:
		- Adenoviridae    - Togaviridae
		- Anelloviridae
		- Arenaviridae
		- Flaviviridae
		- Paramyxoviridae
		- Parvoviridae
		- Picornaviridae
		- Polyomaviridae
		- Poxviridae
		- Rhabdoviridae
			
	'sjf' (Skin or Joint Fluids) includes:
		- Flaviviridae
		- Herpesviridae
		- Papillomaviridae
		- Paramyxoviridae
		- Parvoviridae
		- Picornaviridae
		- Polyomaviridae
		- Rhabdoviridae
		- Togaviridae
				
	'special' (Special pathogens/Zoonotics/Vectorborne) includes:
		- Arenaviridae    - Paramyxoviridae
		- Arteriviridae   - Peribunyaviridae
		- Bornaviridae    - Phenuiviridae
		- Coronaviridae   - Picornaviridae
		- Filoviridae     - Pneumoviridae
		- Flaviviridae    - Poxviridae
		- Hantaviridae    - Reoviridae
		- Hepeviridae     - Retroviridae
		- Herpesviridae   - Rhabdoviridae
		- Nairoviridae    - Togaviridae
'''
	parser = argparse.ArgumentParser(prog='viralstep syndromicbaits', formatter_class=argparse.RawDescriptionHelpFormatter, description='Filter FASTA sequences by virus family or specific groups of families.', usage=USAGE)

	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input fasta file")
	parser.add_argument('-l', '--list', dest="input_family_file", required=True, help="Input corresponding list of virus families identified for each sequence. Columns should be: AccessionID, TaxID, Baltimore classificaton, Virus family.")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="Output filename")
	parser.add_argument ('-f', '--filter', dest="filter_name", required=True, help="Specify virus family (individual or comma-separated list) or pre-defined group of virus families (See usage)")
	parser.add_argument('-v', '--rejected', dest="not_filtered", action="store_true", required=False, help="Obtain remaining sequences (wash-through) rejected by filter. Output contains 'rejected_' prefix in output name")
	parser.add_argument('-t', '--threads', dest="threads", default=None, required=False, help="Number of threads. Will attempt to assign individual threads for each output file. Default = 1.")
	parser.add_argument('--taxid', dest="receive_taxid", action="store_true", required=False, help="Obtain filtered taxID list as additional output. Output contains 'taxID_' prefix in output name")
	parser.add_argument('--lineage', dest="receive_lineage", action="store_true", required=False, help="Obtain filtered database file, derived from input database file, as additional output. Output contains 'family_' prefix in output name")
	parser.add_argument('--log', dest="output_log", action="store_true", required=False, help="Obtain log file output")
	return parser

def load_reference(args):
	accID2taxID = defaultdict(list) # taxID: [accIDs]
	taxID2family = {} # normal dictionary object since one taxID = one family
	family2baltimore = {} # family: baltimore classification (i.e. dsDNA, +ssRNA, etc.)
	family_list = args.input_family_file
	
	# family_list = accID, taxID, baltimore_classification, virus_family
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
	# if ".fasta" not in output: # ensure final output is a .fasta file
		# output = output + ".fasta"
	# else:
		# pass
	with open(output, "w+") as out_tax:
		while True:	
			line_taxID = queue_taxID.get()
			if stop in line_taxID:
				return
			else:
				out_tax.write(line_taxID)

def writer_lineage(output, queue_lineage, stop, args):
	# if ".fasta" not in output: # ensure final output is a .fasta file
		# output = output + ".fasta"
	# else:
		# pass
	with open(output, "w+") as out_lineage:
		while True:	
			line_lineage = queue_lineage.get()
			if stop in line_lineage:
				return
			else:
				out_lineage.write(line_lineage)

def writer_no_filter(output, queue_no_filter, stop, args):
	# if ".fasta" not in output: # ensure final output is a .fasta file
		# output = output + ".fasta"
	# else:
		# pass
	with open(output, "w+") as not_out:
		while True:	
			line_no_filter = queue_no_filter.get()
			if stop in line_no_filter:
				return
			else:
				not_out.write(">" + str(line_no_filter.description) + "\n" + str(line_no_filter.seq).upper() + "\n")

def writer(output, queue_filter, stop, args): # generate output through queues
	# if ".fasta" not in output: # ensure final output is a .fasta file
		# output = output + ".fasta"
	# else:
		# pass
### Output core file (filtered list)		
	with open(output, "w+") as out:
		while True:
			line_filter = queue_filter.get()
			if stop in line_filter:
				return
			else:
				out.write(">" + str(line_filter.description) + "\n" + str(line_filter.seq).upper() + "\n")

def writer_alt(output, queue_filter, queue_no_filter, queue_taxID, queue_lineage, stop, args): # generate output through queues
	# output_file=output
	# if ".fasta" not in output_file: # ensure final output is a .fasta file
		# output_file = output_file + ".fasta"
	# else:
		# pass
### Output core file (filtered list)
### output is a list of file paths to output files [main, lineage, taxID, rejected]
	with open(output[0], "w+") as output:
		if args.receive_taxid is True:
			writer_taxID(output[2], queue_taxID, stop, args)
		else:
			pass
		if args.receive_lineage is True:
			writer_lineage(output[1], queue_lineage, stop, args)
		else:
			pass
		if args.not_filtered is True:
			writer_no_filter(output[3], queue_no_filter, stop, args)
		else:
			pass
		while True:
			line_filter = queue_filter.get()
			if stop in line_filter:
				return
			else:
				output.write(">" + str(line_filter.description) + "\n" + str(line_filter.seq).upper() + "\n")
				
def filter_family(batch, args, threads, output):
### Parameter setup
	accID2taxID, taxID2family, family2baltimore = load_reference(args)
	input_file = batch
	output_file = output[0]
	lineage_file = output[1]
	taxID_file = output[2]
	rejected_file = output[3]
	#output_file=args.output_fasta_file
	# if ".fasta" not in output_file: # ensure final output is a .fasta file
		# output_file = output_file + ".fasta"
	# else:
		# pass

	"""	Parameters
	virus_families = ["Ackermannviridae","Adenoviridae","Adomaviridae","Alloherpesviridae","Alphaflexiviridae","Alphasatellitidae","Alphatetraviridae","Alvernaviridae","Amalgaviridae","Amnoonviridae","Ampullaviridae","Anelloviridae","Arenaviridae",
		"Arteriviridae","Ascoviridae","Asfarviridae","Aspiviridae","Astroviridae","Autolykiviridae","Bacilladnaviridae","Baculoviridae","Barnaviridae","Benyviridae","Betaflexiviridae",
		"Bicaudaviridae","Bidnaviridae","Birnaviridae","Bornaviridae","Bromoviridae","Caliciviridae","Carmotetraviridae","Caulimoviridae","Chrysoviridae","Chuviridae","Circoviridae","Clavaviridae","Closteroviridae","Coronaviridae",
		"Corticoviridae","Cruciviridae","Cruliviridae","Cystoviridae","Deltaflexiviridae","Dicistroviridae","Endornaviridae","Euroniviridae","Filoviridae","Fimoviridae","Flaviviridae","Flexiviridae","Fusariviridae","Fuselloviridae","Gammaflexiviridae",
		"Geminiviridae","Genomoviridae","Globuloviridae","Guttaviridae","Hantaviridae","Hepadnaviridae","Hepeviridae","Herpesviridae","Hypoviridae","Hytrosaviridae","Iflaviridae","Inoviridae","Iridoviridae","Lavidaviridae",
		"Leviviridae","Lipothrixviridae","Luteoviridae","Malacoherpesviridae","Marnaviridae","Marseilleviridae","Medioniviridae","Megabirnaviridae","Mesoniviridae","Metaviridae","Microviridae","Mimiviridae","Mononiviridae","Mymonaviridae","Myoviridae",
		"Mypoviridae","Nairoviridae","Nanoviridae","Narnaviridae","Nimaviridae","Nodaviridae","Nudiviridae","Nyamiviridae","Orthomyxoviridae","Papillomaviridae","Paramyxoviridae","Partitiviridae",
		"Parvoviridae","Peribunyaviridae","Permutotetraviridae","Phasmaviridae","Phenuiviridae","Phycodnaviridae","Picobirnaviridae","Picornaviridae","Pithoviridae","Plasmaviridae",
		"Pleolipoviridae","Pneumoviridae","Podoviridae","Polycipiviridae","Polydnaviridae","Polyomaviridae","Portogloboviridae","Potyviridae","Poxviridae","Qinviridae","Quadriviridae","Reoviridae","Retroviridae","Rhabdoviridae","Roniviridae",
		"Rudiviridae","Sarthroviridae","Secoviridae","Siphoviridae","Smacoviridae","Solemoviridae","Solinviviridae","Sphaerolipoviridae","Sphaeromadae","Spiraviridae","Sunviridae",
		"Tectiviridae","Tobaniviridae","Togaviridae","Tolecusatellitidae","Tombusviridae","Totiviridae","Tristromaviridae","Turriviridae","Tymoviridae","Unclassified","Virgaviridae","Wupedeviridae","Yueviridae"]
	
	human_families = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Astroviridae", "Bornaviridae", "Caliciviridae", "Circoviridae", "Coronaviridae",
		"Filoviridae", "Flaviviridae", "Hantaviridae", "Hepadnaviridae", "Hepeviridae", "Herpesviridae", "Nairoviridae", "Orthomyxoviridae", "Papillomaviridae",
		"Paramyxoviridae", "Partitiviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picobirnaviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae",
		"Poxviridae", "Reoviridae", "Retroviridae", "Rhabdoviridae", "Smacoviridae", "Tobaniviridae", "Togaviridae"]
	
	gastro = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Astroviridae", "Caliciviridae", "Circoviridae", "Coronaviridae", "Filoviridae", "Hantaviridae",
		"Hepadnaviridae", "Hepeviridae", "Nairoviridae", "Orthomyxoviridae", "Partitiviridae", "Parvoviridae", "Picobirnaviridae", "Picornaviridae", "Polyomaviridae",
		"Reoviridae", "Tobaniviridae"]
	
	blood = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepadnaviridae",
		"Hepeviridae", "Herpesviridae", "Nairoviridae", "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae",
		"Pneumoviridae", "Polyomaviridae", "Poxviridae", "Reoviridae", "Retroviridae", "Rhabdoviridae", "Togaviridae"]
	
	respir = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Herpesviridae",
		"Nairoviridae", "Orthomyxoviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae", "Reoviridae", "Togaviridae"]
	
	csf = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepeviridae", "Herpesviridae",
		"Nairoviridae", "Orthomyxoviridae", "Paramyxoviridae", "Parvoviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae", "Rhabdoviridae", "Togaviridae"]
	
	urine = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Flaviviridae", "Paramyxoviridae", "Parvoviridae", "Picornaviridae", "Polyomaviridae", "Poxviridae", "Rhabdoviridae", "Togaviridae"]
	
	sjf = ["Flaviviridae", "Herpesviridae", "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Picornaviridae", "Polyomaviridae", "Rhabdoviridae", "Togaviridae"] # skin or joint fluids
	
	special = ["Arenaviridae", "Arteriviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepeviridae", "Herpesviridae",
		"Nairoviridae", "Paramyxoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Poxviridae", "Reoviridae", "Retroviridae",
		"Rhabdoviridae", "Togaviridae"] # Special pathogens/Zoonotics/Vectorborne

"""
	
### Input check -f or --filter

	if args.filter_name.lower() in ('gastro', 'gastrointestinal'):
		#filter = gastro
		filter = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Astroviridae", "Caliciviridae", "Circoviridae", "Coronaviridae", "Filoviridae", "Hantaviridae",
		"Hepadnaviridae", "Hepeviridae", "Nairoviridae", "Orthomyxoviridae", "Partitiviridae", "Parvoviridae", "Picobirnaviridae", "Picornaviridae", "Polyomaviridae",
		"Reoviridae", "Tobaniviridae"]

	elif args.filter_name.lower() in 'blood':
		#filter = blood
		filter = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepadnaviridae",
		"Hepeviridae", "Herpesviridae", "Nairoviridae", "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae",
		"Pneumoviridae", "Polyomaviridae", "Poxviridae", "Reoviridae", "Retroviridae", "Rhabdoviridae", "Togaviridae"]

	elif args.filter_name.lower() in ('respir','respiratory'):
		#filter = respir
		filter = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Herpesviridae",
		"Nairoviridae", "Orthomyxoviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae", "Reoviridae", "Togaviridae"]

	elif args.filter_name.lower() in 'human':
		#filter = human_families
		filter = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Astroviridae", "Bornaviridae", "Caliciviridae", "Circoviridae", "Coronaviridae",
		"Filoviridae", "Flaviviridae", "Hantaviridae", "Hepadnaviridae", "Hepeviridae", "Herpesviridae", "Nairoviridae", "Orthomyxoviridae", "Papillomaviridae",
		"Paramyxoviridae", "Partitiviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picobirnaviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae",
		"Poxviridae", "Reoviridae", "Retroviridae", "Rhabdoviridae", "Smacoviridae", "Tobaniviridae", "Togaviridae"]

	elif args.filter_name.lower() in 'csf':
		filter = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepeviridae", "Herpesviridae",
		"Nairoviridae", "Orthomyxoviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae", "Rhabdoviridae", "Togaviridae"]

	elif args.filter_name.lower() in 'urine':
		filter = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Flaviviridae", "Paramyxoviridae", "Parvoviridae", "Picornaviridae", "Polyomaviridae", "Poxviridae", "Rhabdoviridae", "Togaviridae"]

	elif args.filter_name.lower() in 'sjf':
		filter = ["Flaviviridae", "Herpesviridae", "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Picornaviridae", "Polyomaviridae", "Rhabdoviridae", "Togaviridae"]

	elif args.filter_name.lower() in 'special':
		filter = ["Arenaviridae", "Arteriviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepeviridae", "Herpesviridae",
		"Nairoviridae", "Paramyxoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Poxviridae", "Reoviridae", "Retroviridae",
		"Rhabdoviridae", "Togaviridae"]

	else: ### Load master list
		virus_families = ["Ackermannviridae","Adenoviridae","Adomaviridae","Alloherpesviridae","Alphaflexiviridae","Alphasatellitidae","Alphatetraviridae","Alvernaviridae","Amalgaviridae","Amnoonviridae","Ampullaviridae","Anelloviridae","Arenaviridae",
		"Arteriviridae","Ascoviridae","Asfarviridae","Aspiviridae","Astroviridae","Autolykiviridae","Bacilladnaviridae","Baculoviridae","Barnaviridae","Benyviridae","Betaflexiviridae",
		"Bicaudaviridae","Bidnaviridae","Birnaviridae","Bornaviridae","Bromoviridae","Caliciviridae","Carmotetraviridae","Caulimoviridae","Chrysoviridae","Chuviridae","Circoviridae","Clavaviridae","Closteroviridae","Coronaviridae",
		"Corticoviridae","Cruciviridae","Cruliviridae","Cystoviridae","Deltaflexiviridae","Dicistroviridae","Endornaviridae","Euroniviridae","Filoviridae","Fimoviridae","Flaviviridae","Flexiviridae","Fusariviridae","Fuselloviridae","Gammaflexiviridae",
		"Geminiviridae","Genomoviridae","Globuloviridae","Guttaviridae","Hantaviridae","Hepadnaviridae","Hepeviridae","Herpesviridae","Hypoviridae","Hytrosaviridae","Iflaviridae","Inoviridae","Iridoviridae","Lavidaviridae",
		"Leviviridae","Lipothrixviridae","Luteoviridae","Malacoherpesviridae","Marnaviridae","Marseilleviridae","Medioniviridae","Megabirnaviridae","Mesoniviridae","Metaviridae","Microviridae","Mimiviridae","Mononiviridae","Mymonaviridae","Myoviridae",
		"Mypoviridae","Nairoviridae","Nanoviridae","Narnaviridae","Nimaviridae","Nodaviridae","Nudiviridae","Nyamiviridae","Orthomyxoviridae","Papillomaviridae","Paramyxoviridae","Partitiviridae",
		"Parvoviridae","Peribunyaviridae","Permutotetraviridae","Phasmaviridae","Phenuiviridae","Phycodnaviridae","Picobirnaviridae","Picornaviridae","Pithoviridae","Plasmaviridae",
		"Pleolipoviridae","Pneumoviridae","Podoviridae","Polycipiviridae","Polydnaviridae","Polyomaviridae","Portogloboviridae","Potyviridae","Poxviridae","Qinviridae","Quadriviridae","Reoviridae","Retroviridae","Rhabdoviridae","Roniviridae",
		"Rudiviridae","Sarthroviridae","Secoviridae","Siphoviridae","Smacoviridae","Solemoviridae","Solinviviridae","Sphaerolipoviridae","Sphaeromadae","Spiraviridae","Sunviridae",
		"Tectiviridae","Tobaniviridae","Togaviridae","Tolecusatellitidae","Tombusviridae","Totiviridae","Tristromaviridae","Turriviridae","Tymoviridae","Unclassified","Virgaviridae","Wupedeviridae","Yueviridae"]
		
		try:
			filter = [fam.strip("\n") for fam in args.filter_name.split(",")]
			for fam in filter:
				if fam.lower() not in [x.lower() for x in virus_families]:
					#filter = [args.filter_name]
					print("Invalid virus family: %s" %(fam))
					filter.remove(fam)
				else:
					pass # no change to list
		except: # test error exceptions
			exit("Something went wrong! Invalid virus families specified!")


		
### Obtain and filter sequences
### Check accessionID and search in accID2taxID for corresponding taxID
### Use obtained taxID to search key in taxID2family
### assign family as search_term
### compare and if it matches, take seq_record
### Output directly to file
	

### Single thread file setup 

	if threads < 2:
		output = open(output_file, "w+")
		if args.receive_taxid is True:
			out_tax = open(taxID_file, "w+")
		else:
			pass
		if args.receive_lineage is True:
			out_lineage = open(lineage_file, "w+")
		else:
			pass
		if args.not_filtered is True:
			not_out = open(rejected_file, "w+")
		else:
			pass
	
	for seq_record in input_file:
		for x,y in accID2taxID.items():
			# isolates AccessionID from any sequence header (including post BaitsTools)
			# Interpreting bait header
			if seq_record.description.split(" ")[0].rsplit("_",1)[0] in y: 
				search_term = str(taxID2family[x]).split("\n")[0] # virus family
			# Interpreting regular accession ID
			# Opposite true if db file contains bait headers
			elif seq_record.description.split(" ")[0] in y: 
				search_term = str(taxID2family[x]).split("\n")[0] # virus family
			else:
				continue

			baltimore_term = family2baltimore[search_term]

			print(search_term) # multiprocessing skews output and syncing it slows it down

			if search_term.lower() in [x.lower() for x in filter]:
				if threads < 2:
					output.write(">" + seq_record.description + "\n" + str(seq_record.seq).upper() + "\n")
				else:
					queue_filter.put(seq_record)

### Optional output
				if args.receive_taxid is True:
					if threads < 2:
						out_tax.write(x + "\n")
					else:
						queue_taxID.put(x + "\n")

				if args.receive_lineage is True:
					if threads < 2:
						#out_lineage.write(x + "," + search_term + "\n")
						out_lineage.write(str(seq_record.id) + "," + x + "," + baltimore_term + "," + search_term + "\n")
					else:
						#queue_lineage.put(x + "," + search_term + "\n")
						queue_lineage.put(str(seq_record.id) + "," + x + "," + baltimore_term + "," + search_term + "\n")
			else:
				if args.not_filtered is True:
					if threads < 2:
						not_out.write(">" + seq_record.description + "\n" + str(seq_record.seq).upper() + "\n")
					else:
						queue_no_filter.put(seq_record)

			# else:
				# pass # push output to log, but continue
	
### Close optional output files
	if threads < 2:
		output.close()
		if args.receive_taxid is True:
			out_tax.close()
		else:
			pass
		if args.receive_lineage is True:
			out_lineage.close()
		else:
			pass
		if args.not_filtered is True:
			not_out.close()
		else:
			pass
				
	
### Check output file
### QC on output FASTA file i.e. check that values are present - no empty file(s)
### Output before and after numbers of sequences

def log_file(args, start, finish):
	with open(args.output_fasta_file + ".log", "w+") as output_log:
		output_log.write("Start: " + start + "\n" + "Finish: " + finish + "\n")
	print("\n" + "Start: " + start)
	print("Finish: " + finish)
	print("Filtered sequences based on virus families. Output file: " + args.output_fasta_file)

	
def run():
	processes = []
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	parser = create_parser()
	args = parser.parse_args()
	input = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(input))
	try:
		out_dir = os.path.abspath(args.output_fasta_file).rsplit("/", 1)[0]
	except IndexError: # same as cwd
		out_dir = os.getcwd()
	output_name = str(os.path.basename(args.output_fasta_file))
	if output_name.startswith(("family_", "taxID_", "rejected_")):
		output_name = output_name.split("_", 1)[1]
	output_path = os.path.join(out_dir, output_name) # guaranteed output
	lineage_path = None
	taxID_path = None
	rejected_path = None
	
	num_seq = count_seq(input)
	expected_output = 1
	if args.receive_taxid is True:
		expected_output = expected_output + 1
		taxID_path = os.path.join(out_dir, "taxID_%s.csv" %(output_name))
	if args.receive_lineage is True:
		expected_output = expected_output + 1
		lineage_path = os.path.join(out_dir,"family_%s.csv" %(output_name))
	if args.not_filtered is True:
		expected_output = expected_output + 1
		rejected_path = os.path.join(out_dir, "rejected_%s" %(output_name))
		
	output_list = [output_path, lineage_path, taxID_path, rejected_path]
	
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
			
### Writer protocol			
			writer_process = Process(target=writer, args = (output_list[0], queue_filter, stop_token, args))
			writer_process.start()
			if args.receive_taxid is True:
				taxID_write = Process(target=writer_taxID, args = (output_list[2], queue_taxID, stop_token, args))
				taxID_write.start()
			if args.receive_lineage is True:
				lineage_write = Process(target=writer_lineage, args = (output_list[1], queue_lineage, stop_token, args))
				lineage_write.start()
			if args.not_filtered is True:
				not_write = Process(target=writer_no_filter, args = (output_list[3], queue_no_filter, stop_token, args))
				not_write.start()

### Process class
			for i in range(len(batch)): # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
				p = Process(target=filter_family, args=(batch[i],args, int(args.threads), output_list))
				processes.append(p)
				p.deamon = True
				p.start()
			for p in processes:
				p.join()			
			
			queue_filter.put(stop_token)
			queue_no_filter.put(stop_token)
			queue_taxID.put(stop_token)
			queue_lineage.put(stop_token)
			
			writer_process.join()
			if args.receive_taxid is True:
				taxID_write.join()
			else:
				pass
			if args.receive_lineage is True:
				lineage_write.join()
				#writer_lineage(output_file, queue_lineage, stop, args)
			else:
				pass
			if args.not_filtered is True:
				not_write.join()
				#writer_no_filter(output_file, queue_no_filter, stop, args)
			else:
				pass
		
		except ValueError: # just give up
			print("Not enough processes for simultaneous writing of output files. Running single process!")
			batch = split_fasta(input, num_seq, int(1))
			filter_family(batch[0], args, int(1), output_list)
				
		except Exception: # one writer process - arises when 1 <= num_seq <= expected_output (too few sequences relative to threads or vice-versa)
			try:
				if int(num_seq-1) <= 0:
					raise ValueError("Only one sequence (and therefore) thread available.")
					#else:
					#	batch = split_fasta(input, num_seq, int(num_seq)-1)
					#	print("Split (uneven 2) input FASTA across %s thread(s), with 1 thread dedicated to writing." % (int(num_seq)-1))
				else:
					args.threads = int(num_seq) 
					batch = split_fasta(input, num_seq, int(args.threads)-1)
					#batch = split_fasta(input, num_seq, int(args.threads)-1) # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
					print("Split input FASTA across %s thread(s), based on number of sequences, with 1 thread dedicated to writing." % (int(args.threads)-1))
				
				writer_process = Process(target=writer_alt, args = (output_list, queue_filter, queue_no_filter, queue_taxID, queue_lineage, stop_token, args))
				writer_process.start()
			
				for i in range(len(batch)): # int(args.threads)-1 (or expected_output) if we want to process with one less thread, save one for writer process
					p = Process(target=filter_family, args=(batch[i],args, int(args.threads), output_list))
					processes.append(p)
					p.deamon = True
					p.start()
				for p in processes:
					p.join()			
			
				queue_filter.put(stop_token)
				queue_no_filter.put(stop_token)
				queue_taxID.put(stop_token)
				queue_lineage.put(stop_token)
			
				writer_process.join()
				
			except ValueError: # just give up
				print("Too few sequences. Running single process!")
				batch = split_fasta(input, num_seq, int(1))
				filter_family(batch[0], args, int(1), output_list)
	
	else: # don't specify thread number
		print("No thread number provided. Running single process!")
		batch = split_fasta(input, num_seq, int(1))
		filter_family(batch[0], args, int(1), output_list)
	
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	
	if args.output_log is True:
		log_file(args, start, finish)
	else:
		print("\n" + "Start: " + start)
		print("Finish: " + finish)
		print("Filtered sequences based on virus families. Output file: " + args.output_fasta_file)

if __name__ == '__main__':
	queue_filter = Queue() # confirm - takes separate process (would explain +1 thread)
	queue_no_filter = Queue()
	queue_taxID = Queue()
	queue_lineage = Queue()
	# manager = Manager()
	# queue_filter = manager.Queue() # confirm - takes separate process (would explain +1 thread)
	# queue_no_filter = manager.Queue()
	# queue_taxID = manager.Queue()
	# queue_lineage = manager.Queue()
	stop_token = "Stop!"
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	try:
		directory = script_path.rsplit("/", 1)[0]
	except IndexError:
		directory = os.getcwd() # script in current working directory
	run()
	print("\n" + "Exited at: " + str(datetime.datetime.utcnow()).split('.')[0])
	exit()