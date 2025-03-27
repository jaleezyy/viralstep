#!/usr/bin/env python3

### Goal of this script is to create a two-pass filter/cluster system
### Use annotation (i.e., Prokka) to determine framework for a given virus family (this using a subset of sequences) and then follow up with mmseq or sequence-based clustering methods
### The point is to determine the expected gene families for a given virus family and then allow sub-sequences to cluster to one of those ORFs
### Example: coronaviruses annotated to ORFs such as RdRP or N or E. Within this initial framework, there will be some sequences in those clusters from the annotated cluster. We can then use MMSEQ using those initial clusters. If score is higher than a threshold, then we cluster (we can check all ORFs and assign the highest score)
### Sequence-based methods alone succumb to high sequence variation leading to high cluster counts with few sequences
### Annotated-based methods alone are too stringent to what we already know and can lead to blind spots because sequence space variation is smaller or not allowed to fit within the expected sequence space
### Annotated + sequence clusters provide a framework but allow sequences similar but may not be annotated as such to group together (i.e., recombination between some other virus family) 

import os, sys, subprocess
from shutil import rmtree, move
from collections import defaultdict
import argparse
from Bio import SeqIO
import graphviz
#from enum import Enum
import random_sequences # custom script
#from random_sequences import random_sample # takes a list of seqs and return subset list

# class SeqContent(Enum):
	# NUCL = 'nucl'
	# PROT = 'prot'

class Cluster:
	"""
	Designate a single cluster representing one ORF. Can add or remove sequences from cluster. No duplicates!
	
	Cluster object:
		ORF Name: Identity of the sequence cluster (i.e., an individual ORF)
		ORF Origin: Taxonomic origin of the ORF (i.e., virus family)
		Sequence: Sequence converted to all upper-case
		
	Expected functions:
		- add sequences
		- remove sequences
		- count number of sequneces within the cluster
	
	"""
	def __init__(self, orf, start, header=None, acc='ref', vf='framework'):
		"""
		For a cluster, define ORF name with starting string sequence or list object of sequences.
		
		Assume header is in the format of <acc>_<virus_family>.
		
		Cluster object is a dict of {acc_vf: set(seqs)}
		"""
		self.orf = orf.lower()
		self.clust = defaultdict(set) 
		# total (unique) sequence cluster
		# accID_vf: set(seqs)
		self.nodes = defaultdict(set)
		# store node/edge connections
		# node (seq_name) -> set(edge/seq_name)

		if header is not None and isinstance(header, list):
			acc = [h.rsplit("_",1)[0].upper() for h in header]
			vf = [v.rsplit("_",1)[1].title() for v in header]
			assert len(acc) == len(vf)
		elif header is not None:
			acc = header.rsplit("_",1)[0].upper()
			vf = header.rsplit("_",1)[1].title()
			assert len(list(acc)) == len(list(vf))
		else:
			pass # keep defaults or what is provided
		if isinstance(start, list):
			seqs = [s.upper() for s in start]
		else:
			seqs = [start.upper()]
		
		if isinstance(acc, list) and len(seqs) == len(acc):
			for a,v,s in zip(acc, vf, seqs):
				self.clust[f"{a}_{v}"].add(s)
		else:
			for s in seqs:
				self.clust[f"{acc}_{vf}"].add(s)
	
	def add_sequence(self, seq, header=None, acc="noAcc", vf="noFamily"):
		"""
		Add a single sequence into existing cluster object
		"""
		seq = str(seq).upper()
		if header is not None:
			try:
				assert len(header.rsplit("_",1)) == 2
			except AssertionError:
				print(header)
				sys.exit(1)
			acc = str(header).rsplit("_",1)[0].upper()
			vf = str(header).rsplit("_",1)[1].title()
		self.clust[f"{acc}_{vf}"].add(seq)
	
	def remove_sequence(self, seq, header=None):
		"""
		Remove a sequence from existing cluster object
		"""
		seq = str(seq).upper()
		# if header provided, we're removing the whole thing
		if header is not None:
			assert len(header.rsplit("_",1)) == 2
			acc = str(header).rsplit("_",1)[0].upper()
			vf = str(header).rsplit("_",1)[1].title()
			self.clust.pop(f"{acc}_{vf}", "Sequences not found!")
		else:
			# find sequence(s) first
			# return {acc}_{vf} where sequence is in the corresponding set
			to_remove = [s for s in self.clust if f"{seq}" in self.clust[s]]
			if len(to_remove) > 1: print("WARNING: Multiple instances of the sequence found in cluster!") 
			# convert to proper WARNING via logging module
			if len(to_remove) > 0:
				for item in to_remove:
					self.clust[item].remove(seq)
		
	def count_sequences(self):
		count = sum([len(self.clust[item]) for item in self.clust])
		# count = 0
		# for item in self.clust:
			# count = count + sum([len(self.clust[item]) for item in self.clust])
		return count
	
	def find_sequence(self, seq, sequence=False):
		"""
		Find sequence within cluster. Returns True, if found
		
		'sequence' indicates whether input 'seq' is a sequence string
		"""
		seq_match = False
		
		if sequence:
			seq = str(seq).upper()
			search = [s for s in self.clust if f"{seq}" in self.clust[s]]
			if len(search) > 0:
				seq_match = True
		else:
			# using name only
			search = [ref for ref in self.clust if f"{seq}" in ref]
			if len(search) > 0:
				seq_match = True
		
		return seq_match
	
	def by_accession(self):
		"""
		Return a form of self.clust where it's a dict of {acc: {seqs}}
		"""
		by_acc = defaultdict(set)
		for item in self.clust:
			acc = str(item).split("_")[0]
			if acc not in by_acc:
				by_acc[acc] = self.clust[item] # initialize set
			else:
				by_acc[acc].update(self.clust[item])
		
		return by_acc
		
	def by_vf(self):
		"""
		Return a form of self.clust where it's a dict of {vf: {seqs}}
		"""
		by_vf = defaultdict(set)
		for item in self.clust:
			vf = str(item).rsplit("_",1)[1].title()
			if vf not in by_vf:
				by_vf[vf] = self.clust[item] # initialize set
			else:
				by_vf[vf].update(self.clust[item])
			
		return by_vf
		
	def by_seq(self):
		"""
		Legacy format where self.clust was just a set of sequences
		"""
		seqs = []
		for item in self.clust:
			seqs.extend([s for s in self.clust[item]])
			
		return set(seqs)
		
	def generate_pairs(self, clust=None):
		"""
		Return list of tuples with each acc and corresponding sequence. 
		
		If alternate clust provided, check if defaultdict (default, by_acc, or by_vf) or set (by_seq) and adjust output accordingly
		"""
		seq_pairs = []
		if clust is None:
			clust = self.clust
		
		if isinstance(clust, defaultdict):
			for id in [i for i in clust]:
				seq_pairs.extend([(id, v) for v in clust[id]])
		elif isinstance(clust, set):
			seq_pairs.extend([(self.orf, v) for v in clust])
			
		return seq_pairs
	
	def print_cluster(self):
		"""
		Use for debugging
		"""
		# sequences
		print(f"For {self.orf}, the following sequences are found:")
		for item in self.clust:
			prev = item[2][0:6]
			print(f">{item[0]}_{item[1]}\n{prev}...\n")
			
		# nodes
	
	def add_node(self, node, edge, label=None):
		"""
		Ensure ref (i.e., Node) is found within cluster by name and then note seq_name as an edge connection
		"""
		ref_name = node
		seq_name = edge
		# check if ref is valid
		# if so, add node
		if self.find_sequence(node):
			self.nodes[node].add(edge)
			
	def remove_node(self, node, edge=None):
		"""
		Check if ref_name (node) in nodes. If seq_name (edge) provided, remove seq, else remove entire node.
		"""
		ref_name = node
		seq_name = edge
		if ref_name in self.nodes:
			if seq_name is not None:
				self.nodes[ref_name].remove(seq_name)
			else:
				self.nodes.pop(ref_name)
				
	def generate_nodes(self):
		"""
		Return list of tuples with (node, edge) pairs.
		"""
		node_pairs = []
		node_terms = {}
		
		# create tuples of (node, edge)
		# this will be used to assign connections
		for node in [i for i in self.nodes]:
			node_pairs.extend([(node, e, None) for e in self.nodes[node]])
		
		# assign each term to a header i.e., TERM1, TERM2, etc
		# will be used to label connections
		counter = 1
		for pair in node_pairs:
			if pair[0] not in node_terms:
				node_terms[pair[0]] = f"TERM{counter}"
				counter+=1
			if pair [1] not in node_terms:
				node_terms[pair[1]] = f"TERM{counter}"
				counter+=1
		
		return node_pairs, node_terms
	
class Framework:
	"""
	Reference framework (ideally for a single virus family). Will contain a series of Cluster objects. Framework should take input FASTA and produce cluster objects.

	1) Random subset will undergo annotation
	2) Parsing annotation will produce Cluster objects with name: set(sequence)
	3) Can pull clustered sequences and produce necessary alignment dbs
	4) Based on match, add sequence to cluster

	Framework --> orf1 = Cluster(orf1, start) --> orf1.add_sequence(seq)
			  --> orf2 = Cluster(orf2, start) --> orf2.add_sequence(seq)
			  --> orf3 = Cluster(orf3, start) --> orf3.add_sequence(seq)
	
	Expected functions:
		- add sequence (upper-level function to call Cluster.add_sequence())
		- remove sequnece (upper-level function to call Cluster.remove_sequence())
		- Visually represent the existing framework with clusters and sequnces within said clusters (likely by name only)
	
	"""
	def __init__(self, ref=None, content='nucl'):
		self.clusters = {} # ORF: Cluster(ORF, seqs, header)
		self.taxonomy = defaultdict(dict) # ORF: {header: seqs} 
		self.ref = ref # starting virus family to build clusters for
		self.content = content
		assert self.content in ('nucl', 'prot')
		
	def add_cluster_object(self, object, header=None):
		"""
		Add one Cluster object into the Framework, adding only the sequences, if cluster exists. Allows for multiple sequences to be input at once.
		"""
		orf_name = object.orf.lower()
		if orf_name not in self.clusters:
			self.clusters[orf_name] = object
		else: 
			for seq in object.clust:
				self.clusters[orf_name].add_sequence(seq, header=header)
	
	def remove_cluster_object(self, name):
		"""
		Remove an entire Cluster object from the Framework
		"""
		self.clusters.pop(name, "Cluster not found!")
		
	def add_cluster_sequence(self, name, seq, header=None, acc=None, vf=None):
		"""
		Add one sequence to existing framework, creating the cluster, if needed. Fine-tune addition of a single sequence. 
		"""
		orf_name = name.lower()
		if acc is not None and vf is not None:
			header = f"{acc.upper()}_{vf.lower()}"
		if orf_name not in self.clusters:
			self.clusters[orf_name] = Cluster(orf_name, seq, header=header)
		else:
			self.clusters[orf_name].add_sequence(seq, header=header)
		
	def remove_cluster_sequence(self, cluster_name, seq, header=None):
		"""
		Call a cluster within a framework and remove a sequence to Cluster object
		"""
		orf_name = cluster_name.lower()
		self.clusters[orf_name].remove_sequence(seq, header)
	
	def pull_cluster_sequences(self, name):
		"""
		Pull all sequences for a given ORF cluster, if found.
		
		Output will be tuples of (acc_id, seq) in a list
		"""
		subset = []
		try:
			orf = self.clusters[name]
			for acc in orf:
				subset.extend([(acc, i) for i in subset.clust[acc]])
		except KeyError:
			subset = [] # reset and exit
			
		return subset
		
	def pull_taxnonomy(self, seq):
		"""
		Given a sequence, use a reverse search to determine the likely taxonomy (i.e., virus family).
		"""
		seq = seq.upper()
		# pull version of each ORF Cluster by virus family
		vf_clusters = [self.clusters[c].by_vf() for c in self.clusters]
		
		search_list = [s for s in vf_clusters if f"{seq}" in vf_clusters[s]]
		
		if len(search_list) == 1 or set(search_list) == 1:
			return list(set(search_list))[0]
		else:
			counter = {}
			for val in search_list:
				if val in counter:
					counter[val] += 1
				else:
					counter[val] = 1
			represent = max(counter, key=counter.get)
			return represent
			
	def add_node_connection(self, cluster_name, node, edge, label=None):
		"""
		Add node to specific cluster within the framework. Will invoke the Cluster.add_node() function.
		
		Returns tuple of (node, edge, label)
		"""
		orf = cluster_name.lower()
		# check if orf exists in sequence set
		# verifies that we have a corresponding cluster prior to registering the node
		if orf not in self.clusters or node not in self.clusters[orf].clust:
			return (None, None, None)
		
		self.clusters[orf].add_node(node, edge, label)
		
		return (node, edge, label)
		
	def remove_node_connection(self, cluster_name, node, edge=None):
		"""
		Remove node from specific cluster within the framework. Will involde the Cluster.remove_node() function.
		"""
		orf = cluster_name.lower()
		# check if orf in node set
		# remove entire node if edge = None
		# if edge specified, the specific node/edge connection will be removed ONLY
		if orf not in self.clusters:
			return (None, None, None)
		
		self.clusters[orf].remove_node(node, edge)
		
		return (node, edge, label)
		
	
	def print_framework(self, output=None, short=True, append=False):
		"""
		Visual representation of existing framework object
		
		Use by_seq() to produce list of sequences found in cluster
		"""
		if (output is not None) and (not append):
			out = open(output, 'w+')
		elif (output is not None) and (append):
			out = open(output, 'a')
		else:
			out = None
		### ADD NODE OUTPUT
		print(f"ORF Clusters:", file=out)
		for clust in self.clusters:
			if short:
				out_list = self.clusters[clust].count_sequences()
				print(f"{clust}: {out_list} sequences", file=out)
			else:
				print(f"{clust}: [{self.clusters[clust].clust}]", file=out)
		
		try:
			out.close()
		except AttributeError: # out = None
			pass
			
		#return ''
		
	def output_framework(self, output_prefix=None, append=False, engine='dot'):
		"""
		Final function to output clusters as individual FASTA file(s). Functions similarly to print_framework; however, separate functions allow for per file generation. 
		
		Prefix required for general filename. Default may be set to 'cluster' for 'ORF_cluster.fasta'. Will overwrite any existing file(s). 
		
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
		allowed_engines = ['dot', 'fdp', 'neato']
		if output_prefix is None:
			print("No output files produced!")
			return None
		else:
			outputs = []
			#output_prefix = str(output_prefix)
			for orf in self.clusters:
				# Collect seqs
				seqs = self.clusters[orf].generate_pairs()
				nodes, labels = self.clusters[orf].generate_nodes()
				
				# Set output files
				out_fasta = f"{orf}_{str(output_prefix)}.fasta"
				output_file = os.path.join(os.getcwd(), out_fasta)
				if (os.path.exists(output_file)) and (not append): 
					os.remove(output_file)
				outputs.append(output_file)
				
				# Sequence-based output
				with open(output_file, 'a') as outclust:
					for seq_pair in seqs:
						outclust.write(f">{seq_pair[0]}\n")
						outclust.write(f"{seq_pair[1]}\n")
						
				# Taxonomic-based output
				# DOT Language per ORF
				out_dot = f"{out_fasta.rsplit('.', 1)[0]}.dot"
				output_dot_file = os.path.join(os.getcwd(), out_dot)
				if (os.path.exists(output_dot_file)) and (not append): 
					os.remove(output_dot_file)
				outputs.append(output_dot_file)
				
				# u007b = "{" and u007d = "}"
				with open(output_dot_file, 'a') as outdot:
					outdot.write(f"digraph {orf.replace('-', '')} " + "\u007b" + "\n") 
					outdot.write(f"rankdir=LR;\n")
					outdot.write(f"splines=true;\n")
					outdot.write(f"overlap=false;\n")
					outdot.write(f"ranksep=10;\n")
					outdot.write(f"node [shape=rect];\n")
					
					# Output individual connections 
					for pair in nodes:
						if pair[2] is not None:
							outdot.write(f"{labels[pair[0]]} -> {labels[pair[1]]} [label={pair[2]};\n")
						else:
							outdot.write(f"{labels[pair[0]]} -> {labels[pair[1]]};\n")

					# Output TERM labels
					for term in labels:
						outdot.write(f"{labels[term]} [label=\"{term}\"];\n")
					outdot.write("\u007d")
					
				
		
		return outputs
		
	def produce_sketch(self):
		"""
		For each cluster in a given framework, link an associated Mash sketch that will correspond with reference sequences for the cluster.
		"""
		start_dir = os.getcwd()
		# initial start with all cluster directories
		mash_dir = os.path.join(start_dir, "mash_sketch")
		ref_files = []
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
		
		
		# divide into cluster directories with one sequence per FASTA
		for orf in self.clusters:
			fasta_files = []
			#num_seqs = self.clusters[orf].count_sequences()
			#seqs = self.clusters[orf].by_seq()
			seqs = self.clusters[orf].generate_pairs()
			num_seqs = len(seqs)

			# aim is to make mash_dir/spike/*.fasta
			try:
				os.mkdir(orf)
			except FileNotFoundError:
				err = os.path.join(os.getcwd(), 'annotation.tsv')
				print(self.clusters[orf].name)
				#print(self.clusters[orf].clust)
				print([o for o in self.clusters])
				exit(f"Check Prokka intermediate file: {err}!")
			for seq, k in zip(seqs, range(1,num_seqs+1)):
				indiv_fasta = os.path.join(orf, f"{seq[0]}_{k}.fasta")
				fasta_files.append(indiv_fasta)
				with open (indiv_fasta, 'w+') as out:
					out.write(f">{seq[0]}\n")
					out.write(f"{seq[1]}")
			
			### TODO: Consider k-mer size for exons/genes, not genomes
			### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7538862/
			### Optimal k-mer length for alignment-free phylogenetics
			mash_ref = os.path.join(start_dir, f"{orf}_reference")
			if self.content == 'prot':
				# default 9
				params = ['mash', 'sketch', '-p', str(threads), '-k', str(kmer_size), '-n', '-o', mash_ref] + fasta_files
			else:
				# default 21; custom 11
				params = ['mash', 'sketch', '-p', str(threads), '-k', str(kmer_size), '-o', mash_ref] + fasta_files
			subprocess.run(params, check=True)
			ref_files.append(f"{mash_ref}.msh")

		os.chdir(start_dir) 
		# return script to outside mash_sketch and plan deletion of mash_dir
		tmp_files.append(mash_dir)
		
		return ref_files

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

def subset_seq(fasta, fasta_name, tmp, size):
	"""
	Given input FASTA and size of the subset, run the random seq script and produce two outputs: the subset and the rest. Temporarily producing DB file(s). 
	"""
	script = random_sequences.__file__
	
	try:
		os.link(fasta, os.path.join(tmp, os.path.basename(fasta)))
	except FileExistsError:
		pass # only applies if splitting between virus families
	os.chdir(tmp)
	out_file = f"subset_{fasta_name}"
	the_rest = f"rejected_subset_{fasta_name}"
	
	params = ['python3', script, '-i', os.path.basename(fasta), '-o', out_file, '-s', str(size), '-v', '--generatedb']
	subprocess.run(params)
	
	os.chdir(start_dir)
	
	return out_file, the_rest

def annotate_sequence(orf_fasta, training_file, content='nucl'):
	"""
	Given a BioPython Seq object, extract the necessary header information and sequence.
	
	Use Prokka to determine taxonomic identity of the sequence (i.e., ORF) and return updated seq_record
	
	prokka --outdir test --prefix annot_test --kingdom "Viruses" --gcode '11' --prodigaltf subset_coronaviridae.fasta.train subset_coronaviridae.fasta
	
	Output of .faa (amino acid) and .ffn (nucleotide) will require header trimming:
	- remove >ID
	- Look for either " protein" or " polyprotein" and take all terms before it
	- If neither found, take whatever is there
	- if len(word) > 1: check for redundant terms like "small" (temp)
	
	Return a dict of {ORF name: Seq}
	"""
	annot_seqs = defaultdict(list)
	protein_terms = ["protein", "polyprotein", "glycoprotein"]
	redundant_terms = ["small"]
	
	# run Prokka on ORF FASTA in tmp_dir
	#training_file = str(orf_fasta).split("_", 1)[1] + ".train"
	subprocess.run(['prokka', '--outdir', 'prokka', '--prefix', 'annotation', '--kingdom', 'Viruses', '--gcode', '11', '--prodigaltf', training_file, orf_fasta])
	if content == 'nucl':
		process = os.path.join('prokka', 'annotation.ffn')
		if check_files(process, os.getcwd()):
			if os.path.getsize(process) == 0: # if empty, rerun without training file
				subprocess.run(['prokka', '--outdir', 'prokka', '--prefix', 'annotation', '--kingdom', 'Viruses', '--gcode', '11', '--force', orf_fasta])
			if os.path.getsize(process) == 0: # after re-running
				exit("No annotated sequences found! Check input and try again!")
			else:
				pass
		else:
			exit("Something went wrong with the annotation! Missing file!")
	elif content == "prot":
		process = os.path.join('prokka', 'annotation.faa')
		if check_files(process, os.getcwd()):
			if os.path.getsize(process) == 0: # if empty, rerun without training file
				subprocess.run(['prokka', '--outdir', 'prokka', '--prefix', 'annotation', '--kingdom', 'Viruses', '--gcode', '11', '--force', orf_fasta])
			if os.path.getsize(process) == 0:
				exit("No annotated sequences found! Check input and try again!")
			else:
				pass
		else:
			print("Something went wrong with the annotation! Missing file!")
			sys.exit(1)
	else:
		print("Something went wrong with the annotation!")
		sys.exit(1)
	
	## TODO: consider replacing "(", ")", "," characters 
	for seq_record in SeqIO.parse(process, 'fasta'):
		header = str(seq_record.description).lower().split(" ", 1)[1]
		sequence = str(seq_record.seq).upper() # may just place entire record 
		
		# acc = str(seq_record.description).lower().split(" ", 1)[0]
		# try:
			# vf = accID2family[acc]
		# except (ValueError, KeyError, TypeError):
			# vf = None
		
		if header.split(" ")[0].lower().startswith("protein"): 
			# exceptions like Protein 3a --> protein-3a
			header = "-".join(i.lower().strip().replace("/", "-") for i in header.split(" "))
		elif len(header.split(" ")) > 1:
			try:
				pos = [header.split(" ").index(j.lower().strip()) for i,j in enumerate(header.split(" ")) if "protein" in j.lower()][0]
				header = "-".join(i.lower().strip().replace("/", "-") for i in header.split(" ")[:pos]) # extract required name
				if header == '': # while terms like Nucleoprotein work, Polyprotein __ seems to error out due to incorrect position
					raise IndexError
			except IndexError:
				header = str(seq_record.description).lower().split(" ", 1)[1] # re-acquire header
				header = "-".join(i.lower().strip().replace("/", "-") for i in header.split(" "))
				# if len(header.strip()) == 0: # if blank still
					# header = "-".join(i.lower().strip().replace("/", "-") for i in str(seq_record.description).split(" ") # take everything
		else:
			pass
			#header = "-".join(i.lower().strip() for i in header.split(" "))
		
		#with open('headers.log', 'a') as out:
		#	print(header, file=out)
		
		annot_seqs[header].append(sequence)
		
	return annot_seqs

def align_sequence(seq, ref, framework, header=None):
	"""
	Given a single sequence as input, along with a list of reference sketch files, run mash dist against each of the cluster reference sketches.
	
	Parse mash dist output to determine lowest distance (equal to an average mash distance per cluster).
	
	Once distances per orf cluster determine, add sequence to most similar Cluster within the given Framework object.
	"""
	scores = {}
	
	# store individual scores
	nodes_per_orf = defaultdict(dict) # {orf: {ref_id: score}} 

	#unclustered = []
	# produce tmp FASTA file given input sequence
	#tmp_fasta = os.path.join(os.getcwd(), 'align.fasta.tmp')
	if header is not None:
		# header is expected to only be an accession ID
		acc_id = header
		try:
			vf = str(accID2family[acc_id]).title()
		except (NameError, KeyError):
			vf = 'sequence'
	else:
		acc_id = 'unknown'
		vf = 'sequence'
	tmp_fasta = os.path.join(os.getcwd(), f"{acc_id}_{vf}.tmp")
	with open(tmp_fasta, 'w+') as out_tmp:
		out_tmp.write(f">{acc_id}_{vf}\n")
		out_tmp.write(f"{seq}") # single sequence input
		
	for file in ref:
		#print(file)
		if check_files(file):
			orf = os.path.basename(file).rsplit(".",1)[0].rsplit("_",1)[0]
			if orf not in scores:
				scores[orf] = 0
			
			# calculate mash distance for given sequence
			print(f"Aligning sequence against {orf}!")
			orf_dist = os.path.join(os.getcwd(), f"{orf}.dist.out")
			
			params = ['mash', 'dist', '-p', str(threads), file, tmp_fasta]
			
			with open (orf_dist, 'w+') as out:
				subprocess.run(params, stdout=out)
			
			
			# score avg mash distance
			with open(orf_dist) as indist:
				line_counter = float(0)
				for line in indist:
					col = line.split("\t")
					assert len(col) == 5
					# expected ref_id = {orf}/{acc_vf}_k.fasta
					ref_id = col[0].rsplit("/",1)[1].rsplit("_",1)[0]
					# expected query_id = {path}/{acc}_{vf}.tmp
					query_id = col[1].rsplit("/",1)[1].rsplit(".",1)[0]
					mash_dist = float(col[2])
					e_val = float(col[3])
					match_hash = col[4]
					
					scores[orf] = scores[orf] + mash_dist
					line_counter+=1
					
					# store ref_ids per ORF
					# once ORF selected, node created linking seq to refs
					nodes_per_orf[orf][ref_id] = mash_dist
					
				scores[orf] = scores[orf] / line_counter
	
	# clean up before end of function
	os.remove(tmp_fasta)
	
	# determine lowest score and add to corresponding cluster
	# ties will not be added at this stage
	print(scores)
	if len(set(scores.values())) == 1:
		print("Unclustered")
		return framework, seq
	else:
		closest_orf = min(scores, key=scores.get)
		print(f"Assigning to {closest_orf}")
		# add sequence to cluster
		framework.add_cluster_sequence(closest_orf, seq, acc=acc_id, vf=vf)
	# add node last to ensure ref_seq check
	ref_to_add = [r for r in nodes_per_orf[closest_orf] if nodes_per_orf[closest_orf][r] < 1]
	for ref in ref_to_add:
		framework.add_node_connection(closest_orf, ref, f"{acc_id}_{vf}")
	return framework, None

def realign_unclustered(unclustered, framework, output=os.getcwd()):
	attempt = 0
	while len(unclustered) > 0:
		attempt+=1
		count_unclustered = len(unclustered)
		print(f"Unclustered sequences found. Updating alignment, attempt {attempt}")
		reference_clusters = framework.produce_sketch()
		
		for seq_record in unclustered:
			seq_to_align = str(seq_record.seq).upper()
			#acc_to_align = str(seq_record.description).split(" ")[0].replace(">", "").rsplit("_", 1)[0]
			acc_to_align = str(seq_record.id).rsplit("_", 1)[0]
			status = align_sequence(seq_to_align, reference_clusters, framework, header=acc_to_align)
			if status is not None:
				pass # stays in unclustered
			else:
				unclustered.remove(seq_record)

		if len(unclustered) == count_unclustered: # no further reduction
			print(f"\n{count_unclustered} unclustered sequence(s) remain after {attempt} attempts! See unclustered.fasta!\n")
			with open(os.path.join(output, "unclustered.fasta"), 'w+') as out:
				for seq_record in unclustered:
					acc_id = str(seq_record.id)
					seq = str(seq_record.seq).upper()
					try:
						vf = str(accID2family[acc_id]).title()
					except (KeyError, NameError):
						vf = 'unknown'
					out.write(f">{acc_id}_{vf}\n")
					out.write(f"{seq}\n")
			break
	
	return framework

def family_hash(db):
	accID2family = {}
	with open(db) as f:
		for line in f:
			col = [item.strip() for item in line.split(",")]
			accID2family[col[0]] = str(col[-1]).title() # first item contains accID, last contains virus family
		
	return accID2family

def run(args):
	fasta_name = os.path.basename(args.input_fasta_file)
	training_file = args.training_file
	### Output directory
	try:
		out_dir = os.path.abspath(args.output_directory)
		os.mkdir(out_dir) # default to os.getcwd() if only name
	except FileExistsError:
		if args.overwrite:
			rmtree(out_dir)
			os.mkdir(out_dir)
		else:
			exit(f"Output directory {args.output_directory} found! Change directory name or use '--overwrite' to force overwrite!")
	### tmp dir
	try:
		tmp_dir = f"{out_dir}_tmp"
		os.mkdir(tmp_dir)
		tmp_files.append(tmp_dir)
	except FileExistsError:
		if args.overwrite:
			rmtree(tmp_dir) #shutil
			os.mkdir(tmp_dir)
			if tmp_dir not in tmp_files:
				tmp_files.append(tmp_dir)
		else:
			exit("Temporary file/directory exists! Check output and use '--overwrite', if required")
	
	### current dir: start_dir
	### Given a FASTA file of input genome sequences, if a ViralSTEP database is provided, we're going to modify all the input per virus family with syndromic filter
	### Rejected flag will output the rest
	### subset and build framework with the virus family and throw the rest after initial framework is done
	virus_pairs = {}
	# call global variable for accID2family, if it is to be used
	global accID2family
	if args.input_family_file is not None:
		accID2family = family_hash(args.input_family_file)
		
		print("Running ViralSTEP Syndromic Filter to separate virus family!")
		viralstep = "/workspace/lab/mcarthurlab/nasirja/viralstep/scripts/syndromic_filter.py" 
		### replace viralstep with proper executable or integration within ViralSTEP
		syn_out = []
		
		# determine individual virus families
		families = set()
		for line in open(args.input_family_file):
			families.add(line.split(",")[-1])
		
		# create separate virus family directories within tmp_dir
		for fam in families:
			fam_name = fam.strip().lower()
			fam_tmp = os.path.join(tmp_dir, fam_name)
			os.mkdir(fam_tmp)
			params = [f"{viralstep}", '-i', str(args.input_fasta_file), '-l', str(args.input_family_file), '-o', os.path.join(fam_tmp, f"{fam_name}.fasta"), '-f', f"{fam_name}", '--rejected', '-t', str(threads)]
			subprocess.run(params)
			
			syn_out.append((os.path.join(fam_tmp, f"{fam_name}.fasta"), os.path.join(fam_tmp, f"rejected_{fam_name}.fasta"))) # (virus family, the rest)
				
			for i in syn_out:
				virus_pairs[i[0]] = i[1]
	else:
		virus_pairs[args.input_fasta_file] = None
		accID2family = None

	### current dir: start_dir
	### determine size of the subset (defined as 25% of the input FASTA)
	for file in virus_pairs:
		os.chdir(start_dir) # ensure a reset if looping
		input = file
		fname = os.path.basename(input)

		other = virus_pairs[input]
		if other is not None:
			tmp_dir = os.path.join(os.getcwd(), os.path.dirname(input)) # all work per virus family isolated to its own tmp_dir
			vf = os.path.basename(os.path.dirname(input)) # tmp_dir per virus family
		else:
			tmp_dir = os.getcwd()
			vf = None
			
		### Create Framework object
		virus_family = Framework(ref=vf, content=args.content) # start up new framework
		logs.append(virus_family) # add to log
		
		num_seq = len(list(SeqIO.parse(input, 'fasta')))
		num_subset = int((num_seq//4) + (num_seq % 4 > 0))
		try:
			smallest_seq = min([len(seq_record.seq) for seq_record in SeqIO.parse(input, 'fasta')])
		except ValueError: 
		# will occur when splitting between virus families and input = 0 sequences
		# TODO: Add this outcome to a log file
			rmtree(os.path.dirname(input))
			#del virus_pairs[input]
			continue # continue loop, removing virus family from further analysis
		print(f"Total size of FASTA: {num_seq}")
		print(f"Smallest genome length found: {smallest_seq}")
		print(f"Total size of expectant subset: {num_subset}")
		if smallest_seq < 20000:
			meta = True
		else:
			meta = False
		subset, remaining = subset_seq(input, fname, tmp_dir, num_subset)
	
		### checkpoint, check if expected output is present
		### Set up Prodigal parameters, if passed
		### Produce or copy Prodigal training file
		if check_files([subset, remaining], tmp_dir):
			if not meta:
				if (training_file is not None) and (os.path.exists(training_file)):
					print("Copying existing training file...")
					os.link(training_file, os.path.join(tmp_dir, f"{subset}.train"))
					os.chdir(tmp_dir)
				else:
					training_file = f"{subset}.train"
					train_cmd = ['prodigal', '-i', f"{subset}", '-t', f"{training_file}"]
					# if meta:
						# train_cmd.append(['-p', 'meta']) # meta flag; default single
					os.chdir(tmp_dir)
					print("Creating training file...") 
					subprocess.run(train_cmd)
					if not check_files(training_file):
						print("Missing training file...")
						sys.exit(1)
			else:
				print("Genomes detected are too small for training! Skipping and running Prodigal using Anonymous mode!")
				os.chdir(tmp_dir)
		else:
			print("Something went wrong! Check temporary files/directories and try again!")
			sys.exit(1)

		### current dir: tmp_dir
		### Take subset and run Prodigal followed by Prokka
		pro_cmd = ['prodigal', '-i', str(subset), '-t', f"{subset}.train", '-a', f"aa_{subset}", '-d', f"nucl_{subset}", '-g', '11', '-o', f"{subset}.out"]
		if meta:
			del pro_cmd[3:5] # -t, training file
			pro_cmd.extend(['-p', 'meta'])
		subprocess.run(pro_cmd)
		
		if args.content == "nucl":
			annot_seqs = annotate_sequence(f"nucl_{subset}", f"{subset}.train", args.content)
		elif args.content == "prot":
			annot_seqs = annotate_sequence(f"aa_{subset}", f"{subset}.train", args.content)
		else:
			exit("Something went wrong selecting content! Check arguments for Prodigal or Prokka!")
			#annot_seqs = annotate_sequence(subset)
		

		### current dir: tmp_dir
		### Given annotated sequences, begin constructing the framework
		### For the subset of sequences, create reference mash sketches
		### sequence header (acc_vf) = starting_framework
		for orf in annot_seqs:
			virus_family.add_cluster_object(Cluster(orf, annot_seqs[orf], vf=vf))
			# each input Cluster defaults to ref_{vf} header
		reference_clusters = virus_family.produce_sketch()
		print(f"\nProduced {len(reference_clusters)} reference clusters!\n")
	
		### current dir: tmp_dir
		### Using remaining sequences, run Prodigal using the same training file followed by mash dist per sequence against each reference cluster
		### Determine lowest overall distance score (closest match) and add sequence to cluster
		pro_cmd = ['prodigal', '-i', str(remaining), '-t', f"{subset}.train", '-a', f"aa_{remaining}", '-d', f"nucl_{remaining}", '-g', '11', '-o', f"{remaining}.out"]
		if meta:
			del pro_cmd[3:5] # -t, training file
			pro_cmd.extend(['-p', 'meta'])
		subprocess.run(pro_cmd)
		unclustered = []
		#virus_family.print_framework("before.log")
		
		if args.content == "nucl":
			for seq_record in SeqIO.parse(f"nucl_{remaining}", 'fasta'):
				seq_to_align = str(seq_record.seq).upper()
				# acc_to_align = str(seq_record.description).split(" ")[0].replace(">", "").rsplit("_", 1)[0]
				acc_to_align = str(seq_record.id).rsplit("_", 1)[0]
				virus_family, status = align_sequence(seq_to_align, reference_clusters, virus_family, header=acc_to_align)
				if status is not None:
					#unclustered.append(status)
					assert str(seq_record.seq).upper() == status.upper()
					unclustered.append(seq_record)
		elif args.content == "prot":
			for seq_record in SeqIO.parse(f"aa_{remaining}", 'fasta'):
				seq_to_align = str(seq_record.seq).upper()
				#acc_to_align = str(seq_record.description).split(" ")[0].replace(">", "").rsplit("_", 1)[0]
				acc_to_align = str(seq_record.id).rsplit("_", 1)[0]
				virus_family, status = align_sequence(seq_to_align, reference_clusters, virus_family, header=acc_to_align)
				if status is not None:
					#unclustered.append(status)
					assert str(seq_record.seq).upper() == status.upper()
					unclustered.append(seq_record)
		else:
			exit("Something went wrong selecting content! Check arguments for Prodigal or Mash!")
		
		### OPTIONAL
		### current dir: tmp_dir
		### Using the other sequences (i.e., not belonging to a virus family), run Prodigal to propose ORFs without the training file (we can't assume the same structure)
		### Run align_sequence to use Mash distance alignment against reference clusters
		### TODO?: Split per virus family and run each family individually...
		if other is not None:
			print(f"Attempting to align all other sequences against {os.path.dirname(other)}!")
			reference_clusters = virus_family.produce_sketch()
			oname = os.path.basename(other) # tmp_dir/{oname}
			pro_cmd = ['prodigal', '-i', str(oname), '-a', f"aa_{oname}", '-d', f"nucl_{oname}", '-g', '11', '-o', f"{oname}.out"]
			if meta:
				pro_cmd.extend(['-p', 'meta'])
			subprocess.run(pro_cmd)
			if args.content == "nucl":
				for seq_record in SeqIO.parse(f"nucl_{oname}", 'fasta'):
					seq_to_align = str(seq_record.seq).upper()
					#acc_to_align = str(seq_record.description).split(" ")[0].replace(">", "").rsplit("_", 1)[0]
					acc_to_align = str(seq_record.id).rsplit("_", 1)[0]
					virus_family, status = align_sequence(seq_to_align, reference_clusters, virus_family, header=acc_to_align)
				if status is not None:
					#unclustered.append(status)
					assert str(seq_record.seq).upper() == status.upper()
					unclustered.append(seq_record)
			elif args.content == "prot":
				for seq_record in SeqIO.parse(f"aa_{oname}", 'fasta'):
					seq_to_align = str(seq_record.seq).upper()
					#acc_to_align = str(seq_record.description).split(" ")[0].replace(">", "").rsplit("_", 1)[0]
					acc_to_align = str(seq_record.id).rsplit("_", 1)[0]
					virus_family, status = align_sequence(seq_to_align, reference_clusters, virus_family, header=acc_to_align)
				if status is not None:
					#unclustered.append(status)
					assert str(seq_record.seq).upper() == status.upper()
					unclustered.append(seq_record)
			else:
				exit("Something went wrong selecting content! Check arguments for Prodigal!")
				
		### OPTIONAL
		### Re-create mash sketch and re-run mash dist for any unclustered sequences
		### May be a fixed attempt limit or loop until no further reductions found
		if not args.realign:
			count_unclustered = len(unclustered)
			print(f"{count_unclustered} unclustered sequence(s) remain! Re-run with `--realign-unclustered' to potentially reduce these sequences! See unclustered.fasta!")
			with open(os.path.join(out_dir, "unclustered.fasta"), 'w+') as out:
				for seq_record in unclustered:
					acc_id = str(seq_record.id)
					seq = str(seq_record.seq).upper()
					try:
						vf = str(accID2family[acc_id]).title()
					except (KeyError, NameError):
						vf = 'unknown'
					out.write(f">{acc_id}_{vf}\n")
					out.write(f"{seq}\n")
		else:
			virus_family = realign_unclustered(unclustered, virus_family, out_dir)
		
		### current (starting) dir: tmp_dir
		### With final Framework object containing Cluster objects, produce final output
		### Each ORF will output its own FASTA outside tmp_dir
		### start_dir/{orf}.fasta
		### final dir: start_dir
		os.chdir(out_dir)
		print("Generating FASTA files per cluster generated!")
		if other is None: # no virus family splitting
			output_files = virus_family.output_framework('cluster')
		else:
			os.mkdir(os.path.basename(os.path.dirname(other)))
			if os.path.exists("unclustered.fasta"):
				move("unclustered.fasta", os.path.basename(os.path.dirname(other)))
			os.chdir(os.path.basename(os.path.dirname(other)))
			output_files = virus_family.output_framework('cluster')
		os.chdir(start_dir)
	
	### End of above loop
	### Clean up temp files
	if (args.cleanup) and (len(tmp_files)) > 0:
		print("Cleaning up!")
		for file in tmp_files:
			if os.path.isdir(file): rmtree(file)
			else:
				try:
					os.remove(file)
				except FileNotFoundError:
					pass

	### Standard output
	print(f"\nClusters generated:")
	for i,j in zip(range(1, len(logs)+1), logs):
		print(f"Cluster {i}:")
		j.print_framework()
	
	### Produce any log files (should be found in out_dir)
	if (args.log) and len(logs) > 0:
		for l in logs:
			#l.print_framework() # stdout counts
			l.print_framework(os.path.join(out_dir, "cluster_framework.log"), short=False, append=True) # file output seqs
		print("All clusters and sequences saved in 'cluster_framework.log'!")
	exit()

def test(db=None):
	global accID2family
	if args.input_family_file is not None:
		accID2family = family_hash(args.input_family_file)
	else:
		accID2family = None
	example_orf = "nucleocapsid"
	example_seq = "AAAAGGGGCCCCTTTT"
	example_add_seq = "ATCGATCGATCGATCG"
	example_list = ["AGGGGCGCGCGCGCG", "TAATATTATATA", "GCGAGAGAGGCG"]
	example_list_headers = ['testing_1', 'testing_2', 'testing_3']
	frame_test = Framework()
	cluster = Cluster(example_orf, example_seq, 'testing_single')
	cluster_list = Cluster("nucleocapsid", example_list, example_list_headers)
	print("Checking single sequence input!")
	print(cluster.count_sequences())
	print(cluster.name)
	print(cluster.clust)
	print("Checking list of sequence input!")
	print(cluster_list.count_sequences())
	print(cluster_list.name)
	print(cluster_list.clust)
	frame_test.add_cluster_object(cluster)
	frame_test.add_cluster_object(cluster_list)
	frame_test.add_cluster_sequence("envelope", "GGGCCCGGGCGGC", 'header1_vf')
	frame_test.add_cluster_sequence("nucleocapsid", "AAATATATATAATA", 'header2_vf')
	frame_test.print_framework(short=False) # show all sequences
	frame_test.print_framework(short=True) # default (show only counts per cluster)

def create_parser():
	parser = argparse.ArgumentParser(description='Perform alignment based clustering based on predicted open reading frames (ORFs). If taxonomy metadata provided in ViralSTEP format, the clustering will be broken down per virus family.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file")
	parser.add_argument('-o', '--outdir', dest='output_directory', required=False, default='clusters', help='Output directory containing FASTA files per determined ORF.')
	parser.add_argument('-d', '--database', dest='input_family_file', required=False, default=None, help="ViralSTEP classificaton database file. Will be used to isolate sequences between virus families to produce a framework.")
	parser.add_argument('-c', '--content', dest='content', choices=['nucl', 'prot'], default='nucl', help="Seqeunce content for design. Choose between nucleotide (nucl) or amino acid/protein (prot). Default is 'nucl'.")
	parser.add_argument('-p', '--prodigaltf', dest='training_file', default=None, help="Use existing Prodigal training file. Useful for replicating clusters.")
	parser.add_argument('-r', '--realign-unclustered', dest='realign', action='store_true', help="Update mash sketches and realign unclustered sequences. Will loop until no further reductions are possible. If any unclustered sequences remain, they will be output as a FASTA file.")
	parser.add_argument('-k', '--kmer', dest='kmer', type=int, default=21, help='K-mer size for MASH alignment. Default is 21, matching default parameters for MASH')
	parser.add_argument('--overwrite', dest='overwrite', action='store_true', help="Overwrite any existing intermediate or temporary files/directories")
	parser.add_argument('--cleanup', dest='cleanup', action='store_true', help='Delete temporary directories or files produced during the run.')
	parser.add_argument("--test", dest="test", action='store_true', help="Run test function to ensure basic Cluster and Framework objects are operational. Will exit after running!")
	parser.add_argument('--log', dest='log', action='store_true', help='Produce log file(s) with any Cluster frameworks produced')
	parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads. Default is 1.")
	return parser
	

if __name__ == '__main__':
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	tmp_files = []
	logs = []
	try:
		script_dir = script_path.rsplit("/", 1)[0]
	except IndexError:
		script_dir = os.getcwd() # script in current working directory
	
	#test()
	#print("No errors!")

	parser = create_parser()
	args = parser.parse_args()
	if args.test:
		test(args.input_family_file)
		try:
			print(accID2family)
		except NameError:
			print("Not a global variable!")
		exit()
	global threads
	threads = args.threads
	global kmer_size
	kmer_size = args.kmer
	run(args)
	os.chdir(start_dir)
	exit()