#!/usr/bin/env python3

### Goal of this script is to take a cluster FASTA file and use mash to determine the inter-virus family connections
### Will rely on functions from 'dot_cluster_framework.py'
### Will effectively be a shorter version of the clustering script with the assumption that each sequence does belong regardless of result
### Per sequence, use MASH to align and assign node --> edge connections if mash distance < 1 (remember mash distance of 1 == too divergent)
### Lowest mash distance will be added and mash will re-run
### Connections already made will be ignored and new connections will be formed
### Repeat until all sequences have been processed

### Alternate, possibly faster
### Generate mash sketch with given FASTA
### Run each individual sequence against sketch to score each alignment?
### Generate node --> edge based on mash distance per sequence

import os, sys
import argparse
from shutil import rmtree
from Bio import SeqIO
from collections import defaultdict
import subprocess

class Node:
	"""
	Node object to store and visualize sequence alignment relationships
	"""
	def __init__(self, fasta):
		# establish all nodes upon initialization
		# expected header: {acc_id}_{vf}_{orf}
		self.headers = [str(seq_record.description).rsplit("_", 1)[0] for seq_record in SeqIO.parse(fasta, 'fasta')]
		#self.nodes = defaultdict(set,{ k: set() for k in self.headers })
		self.nodes = defaultdict(set)
		
		# infer ORF name from headers
		self.orf_name = '-'.join(list(set([str(seq_record.description).rsplit("_", 1)[1] for seq_record in SeqIO.parse(fasta, 'fasta')])))
		
	def add_node(self, node, edge, dist=None):
		"""
		Add node --> edge connections
		
		Node must exist as per __init__
		"""
		assert node in self.headers
		
		# if they match, ignore
		if node == edge:
			return None
		
		# if edge has a node, let's check for existing relationships
		# first check is if edge already has a node
		if edge in self.nodes:
			# pull value of matching node
			check = self.nodes[edge] # output should be a set for the given node
			# if reverse edge --> node not in cluster, add it (new connection)
			# if the value for node is not found in the edge equivalent, add new node --> edge
			if node not in [e[0] for e in check]:
				self.nodes[node].add((edge, dist))
			# if edge --> node already exists, but dist differs, add it
			else:
				for e in check:
					if e[0] == node:
						if e[1] != dist:
							self.nodes[node].add((edge, dist))
							break
		else:
			# this is the first instance of this relationship, just add it
			self.nodes[node].add((edge, dist))
		
	def remove_node(self, node, edge=None):
		"""
		Remove an entire node, or a specific connection if edge is provided
		"""
		if edge is not None:
			for e in self.nodes[node]:
				if self.nodes[e][0] is edge:
					self.nodes[node].remove(e)
		else:
			self.nodes.pop(node)

	def generate_nodes(self):
		"""
		Return list of tuples with (node, edge, dist) pairs.
		"""
		# node_pairs = each connection with dist (node, edge, dist)
		node_pairs = []
		# node_terms = each node or edge will be given a numbered label
		node_terms = {}
		
		# create tuples of (node, edge, dist)
		# this will be used to assign connections
		for node in [i for i in self.nodes]:
			node_pairs.extend([(node, e[0], e[1]) for e in self.nodes[node]])
		
		# assign each term to a header i.e., TERM1, TERM2, etc
		# will be used to label connections
		counter = 1
		for pair in node_pairs:
			if pair[0] not in node_terms:
				node_terms[pair[0]] = f"TERM{counter}"
				counter+=1
			if pair[1] not in node_terms:
				node_terms[pair[1]] = f"TERM{counter}"
				counter+=1
				
		return node_pairs, node_terms
		
	def output_nodes(self, filename='output'):
		"""
		Final function to output a DOT file.
		
		Filename required. Will overwrite any existing file(s). 
		
		Node output sample:
		
		digraph orf {
		rankdir=TB;
		node [shape=record];
		TERM1 -> TERM2 [label="label"];
		TERM3 -> TERM4;
		TERM5 -> TERM6;
		TERM7 -> TERM2;
		TERM6 [label="test3_3"];
		TERM1 [label="test1"]; (node)
		TERM2 [label="test1_2"]; (edge)
		TERM7 [label="test4"];
		TERM4 [label="test2_2"];
		TERM5 [label="test3"];
		TERM3 [label="test2"];
		}
		"""
		#allowed_engines = ['dot', 'fdp', 'neato']
		
		for node in self.nodes:
			# Collect node tuples
			nodes, labels = self.generate_nodes()
			
			# Set output files
			if filename.endswith(".dot"):
				out = f"{filename}"
			else:
				out = f"{filename}.dot"
			output_file = os.path.join(os.getcwd(), out)
			if (os.path.exists(output_file)):
				os.remove(output_file)
				
			# Taxonomic-based output
			# DOT Language per ORF
			# u007b = "{" and u007d = "}"
			with open(output_file, 'a') as outdot:
				outdot.write(f"digraph {self.orf_name.replace('-', '')} " + "\u007b" + "\n") 
				outdot.write(f"rankdir=LR;\n")
				outdot.write(f"splines=true;\n")
				outdot.write(f"overlap=false;\n")
				outdot.write(f"ranksep=10;\n")
				outdot.write(f"node [shape=rect];\n")
				
				# Output individual connections 
				for pair in nodes:
					if pair[2] is not None:
						outdot.write(f"{labels[pair[0]]} -> {labels[pair[1]]} [label={pair[2]}];\n")
					else:
						outdot.write(f"{labels[pair[0]]} -> {labels[pair[1]]};\n")

				# Output TERM labels
				for term in labels:
					outdot.write(f"{labels[term]} [label=\"{term}\"];\n")
				outdot.write("\u007d")
				
		return output_file
				

class Align:
	"""
	Object represents a growing cluster that can invoke alignment and return an object of sequence alignment relationships. 
	"""
	def __init__(self, fasta, name=None):
		"""
		For a given fasta, define name for file names
		
		"""
		self.sequences = [] # stored as a list of (acc, seq)
		
		for seq_record in SeqIO.parse(fasta, 'fasta'):
			self.sequences.append((seq_record.description, str(seq_record.seq).upper()))
			
		if name is None:
			self.ref_name = f"{os.path.basename(fasta)}"
		else:
			self.ref_name = name
		
	
	def produce_sketch(self, kmer_size=21):
		"""
		Take all cluster sequences (self.sequences) and generate a mash sketch by splitting and creating individual FASTA files as input.
		"""
		start_dir = os.getcwd()
		# initial start with all cluster directories
		mash_dir = os.path.join(start_dir, "mash_sketch")
		ref_name = self.ref_name

		try:
			os.mkdir(mash_dir)
			#tmp_files.append(mash_dir)
			os.chdir(mash_dir)
		except FileExistsError:
			rmtree(mash_dir)
			os.mkdir(mash_dir)
			#if mash_dir not in tmp_files:
			#	tmp_files.append(mash_dir)
			os.chdir(mash_dir)
			
		# split sequences into mash_dir/{header}_k.fasta
		# append individual files into mash sketch
		fasta_files = []
		num_seqs = len(self.sequences)
		for seq, k in zip(self.sequences, range(1,num_seqs+1)):
			indiv_fasta = os.path.join(mash_dir, f"{seq[0]}_{k}.fasta")
			fasta_files.append(indiv_fasta)
			with open(indiv_fasta, 'w+') as out:
				out.write(f">{seq[0]}\n")
				out.write(f"{seq[1]}")
			
		### TODO: Consider k-mer size for exons/genes, not genomes
		### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7538862/
		### Optimal k-mer length for alignment-free phylogenetics
		mash_ref = os.path.join(start_dir, ref_name)
		params = ['mash', 'sketch', '-p', str(threads), '-k', str(kmer_size), '-n', '-o', mash_ref] + fasta_files
		
		# if content == 'prot':
			# default 9
			# params = ['mash', 'sketch', '-p', str(threads), '-k', str(kmer_size), '-n', '-o', mash_ref] + fasta_files
		# else:
			# default 21; custom 11
			# params = ['mash', 'sketch', '-p', str(threads), '-k', str(kmer_size), '-o', mash_ref] + fasta_files
		
		subprocess.run(params, check=True)

		# return script to outside mash_sketch and plan deletion of mash_dir
		tmp_files.append(mash_dir)
		os.chdir(start_dir) 
		
		return f"{mash_ref}.msh"
		
	
	def align_sequence(self, mash_ref, dot=None):
		"""
		Given a single mash sketch reference file as input, run mash dist using individual sequences from self.sequences.
		
		Parse mash dist output to determine lowest distance (equal to an average mash distance per cluster).
		
		Add nodes if mash distance < 1 for a given sequence and sketch reference.
		
		If dot object provided, add and update nodes, otherwise only output will be a table of all mash output
		
		"""
		mash_output = f"{self.ref_name}.all.dist.out"
		if os.path.exists(mash_output):
			os.remove(mash_output)
			all_dist = open(mash_output, 'w+')
		else:
			all_dist = open(mash_output, 'w+')
			
		# provide headers for output file
		all_dist.write(f"Reference (Node)\tQuery (Edge)\tMash-Distance\tE-Value\tMatching-Hash\n")
		
		
		# produce tmp FASTA file given sequences
		# alternate: refer to mash_dir to pull individual FASTAs
		for seq in self.sequences:
			tmp_fasta = os.path.join(os.getcwd(), f"{seq[0]}.tmp")
			with open(tmp_fasta, 'w+') as out_tmp:
				out_tmp.write(f">{seq[0]}\n")
				out_tmp.write(f"{seq[1]}") # single sequence input
				
			# using provided mash_ref, perform mash dist per sequence
			# store score per reference comparison
			# scores = {} # ref: score
			
			# calculate mash distance for given sequence
			print(f"Aligning {seq[0]} against the cluster!")
			tmp_dist = os.path.join(os.getcwd(), f"{seq[0]}.dist.tmp")
			params = ['mash', 'dist', '-p', str(threads), mash_ref, tmp_fasta]
				
			with open (tmp_dist, 'w+') as out:
				subprocess.run(params, stdout=out)
				
			# score avg mash distance
			with open(tmp_dist) as indist:
				# line_counter = float(0)
				for line in indist:
					col = line.split("\t")
					assert len(col) == 5
					# expected ref_id = {orf}/{acc_id}_{vf}_{orf}_{k}.tmp
					ref_id = col[0].rsplit("/",1)[1].rsplit(".",1)[0].rsplit("_",2)[0]
					# expected query_id = {path}/{acc_id}_{vf}_{orf}.tmp
					query_id = col[1].rsplit("/",1)[1].rsplit(".",1)[0].rsplit("_", 1)[0]
					mash_dist = float(col[2])
					e_val = float(col[3])
					match_hash = col[4]
					
					# output mash line to central output file
					all_dist.write(f"{ref_id}\t{query_id}\t{mash_dist}\t{e_val}\t{match_hash}")
					
					# if dot provided, check mash_dist if < 0
					# update or add to node
					if dot is not None and mash_dist < 1:
						print(f"Adding {ref_id} node with {query_id} edge!")
						dot.add_node(ref_id, query_id, mash_dist)
					
				# scores[orf] = scores[orf] / line_counter
		
			# clean up before end of loop
			os.remove(tmp_fasta)
			os.remove(tmp_dist)
		
		all_dist.close()
		
		return mash_output

def check_files(filename, tmp=None):
	"""
	Use to check if one (string) or more (list) files exists in the expected directory
	"""
	if isinstance(filename, list):
		checklist = {}
		for file in filename:
			if tmp is not None:
				check = os.path.join(tmp, file)
			else:
				check = file # assume current directory i.e., filename
			if os.path.exists(check):
				checklist[file] = True
			else:
				checklist[file] = False
		if False in checklist.values():
			return False
		else:
			return True
	else:
		if tmp is not None:
			check = os.path.join(tmp, filename)
		else:
			check = filename
		if os.path.exists(check):
			return True
		else:
			return False

def cleanup(tmp):
	if len(tmp) < 1:
		return None
	else:
		for item in tmp:
			if os.path.isdir(item):
				rmtree(item)
			else:
				os.remove(item)

def run(args):
	### determine name
	fasta_name = os.path.basename(args.input_fasta_file)
	if args.output_name is None:
		output_name = fasta_name
	else:
		output_name = args.output_name
	
	### invoke align and node objects
	align_seq = Align(args.input_fasta_file, output_name)
	node_seq = Node(args.input_fasta_file)
	
	### Produce mash sketch
	mash_sketch = align_seq.produce_sketch(args.kmer)
	
	### Run mash alignment with supplied Dot object
	final_mash_output = align_seq.align_sequence(mash_sketch, node_seq)
	
	### From updated dot object, produce output
	final_dot_output = node_seq.output_nodes(output_name)
	
	### Assert outputs are found
	assert check_files([final_mash_output, final_dot_output])
	
	### Clean up
	cleanup(tmp_files)

def create_parser():
	parser = argparse.ArgumentParser(description='Perform quick mash alignment per sequential sequence addition to a growing sketch from a cluster FASTA. Generate DOT output')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input cluster FASTA file. Currently works with nucleotide FASTA.")
	parser.add_argument('-o', '--output-name', dest='output_name', default=None, help="Name of ORF cluster (ex. Spike). Default will use the basename of the input fasta file. Will be used in all output files")
	parser.add_argument('-k', '--kmer', default=21, help='K-mer size used for MASH')
	parser.add_argument('-t', '--threads', default=1, help='Number of threads. Default = 1')
	# parser.add_argument('--keep-mash', action='store_true', help="Keep all mash alignment output. Will be stored in a single file in the format: '{name}.all.dist.out', with 'name' referring to the parameter -n or --name. If file already exists, IT WILL BE OVERWRITTEN")
	return parser

if __name__ == '__main__':
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	global threads
	global tmp_files
	try:
		script_dir = script_path.rsplit("/", 1)[0]
	except IndexError:
		script_dir = os.getcwd() # script in current working directory
	
	parser = create_parser()
	args = parser.parse_args()
	threads = args.threads
	tmp_files = []

	run(args)
	sys.exit(0)