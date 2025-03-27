#!/usr/bin/env python

### Goal is to take a FASTA file (or directory of files) and produce a JSON containing ancestral nodes
### Input is presumed to be ORF/gene family sequences
### Each FASTA will undergo MAFFT alignment
### orf_correction.py will eliminate invalid characters first
### trimal will perform automatic trimming --> FastTree
### Run length check to ensure sequences divisible by 3
### orf_correction.py will be re-run to eliminate all stop codons

### length checks will be done using orf_frame_count.py, and run to trim sequences as needed

from Bio import SeqIO
import os, sys, subprocess, shutil, glob
import argparse
import json

def count_char(input, char):
	counter = 0
	for line in open(input):
		counter = counter + line.count(f"{char}")
	
	return counter

def create_parser():
	parser = argparse.ArgumentParser(description='Take FASTA clusters and generate phylogenetic alignments/trees and perform ancestral reconstruction using HyPhy developmental Ancestral Reconstuction software (v0.2).')
	parser.add_argument('-i', '--input', dest="input_files", required=True, help="Input FASTA file or directory of FASTAs. Will auto-detect content.")
	parser.add_argument('-o', '--output', dest="output_dir", default="corrected", help="Output directory containing ancestral nodes")
	parser.add_argument('--keep-intermediate', dest="intermed_keep", action="store_true", required=False, help="Keep intermediate temporary files produced in the run") ### will delete at end of program otherwise
	parser.add_argument('--overwrite', dest='overwrite', action='store_true', help='Overwrite any prior existing intermediate files. If not provided, and intermediates found, they will be used instead of generating new file(s)')
	return parser

def check_scripts(dir):
	"""
	Ensure that required scripts are found in the same directory as this script
	"""
	required_scripts = ['orf_correction.py', 'orf_frame_count.py', 'FitMG94.bf', 'AncestralSequences.bf']
	for script in required_scripts:
		if not os.path.exists(os.path.join(dir, script)):
			exit(f"Required {script} not found!")
			
	### TODO: add check for mafft, trimal and fasttree, HyPhy and corresponding executable files

def mafft(input, name, overwrite=False):
	"""
	Perform gloabal alignment of ORF sequences to later generate a phylogenetic tree
	"""
	output = os.path.join(os.getcwd(), f"{name.rsplit('.', 1)[0]}.align.tmp")
	if os.path.exists(output) and not overwrite:
		print(f"Prior output found! File: {output}")
		return output
	else:
		mafft_cmd = ['mafft', '--maxiterate', '1000', '--globalpair', f"{input}"] 
		
		with open(output, 'w+') as out:
			subprocess.run(mafft_cmd, stdout=out)
	
		assert os.path.exists(output)
		tmp_files.append(output)
	
		return output

def correct_seqs(input, all_stops=False, overwrite=False):
	"""
	Use custom python script to check FASTA file sequences for errors, including the removal of all stop codons and ambiguous characters
	"""
	correction = os.path.join(script_path, 'orf_correction.py')
	output = os.path.join(os.getcwd(), "corrected_"+str(os.path.basename(input))) 
	if os.path.exists(output) and not overwrite:
		print(f"Prior output found! File: {output}")
		return output
	else:
		correction_cmd = [f"{correction}", '-i', f"{input}", '-o', 'corrected', '-c', 'nucl']
		if all_stops:
			correction_cmd.append('--all-stops')
		
		subprocess.run(correction_cmd)
			
		assert os.path.exists(output)
		tmp_files.append(output)
		
		return output

def trim_align(input, overwrite=False):
	"""
	Use trimal to perform automated trimming of a global alignment prior to tree generation
	"""
	output = os.path.join(os.getcwd(), "trimmed_"+str(os.path.basename(input)))
	if os.path.exists(output) and not overwrite:
		print(f"Prior output found! File: {output}")
		return output
	else:
		trim_cmd = ['trimal', '-in', f"{input}", '-out', f"{output}", '-fasta', '-automated1']
		
		subprocess.run(trim_cmd)
	
		assert os.path.exists(output)
		tmp_files.append(output)
	
		return output
		
def correct_length(input, overwrite=False):
	"""
	Use custom python script to ensure the lengths of each sequence in a given file is divisible by 3 (i.e., ensure all codons can be derived)
	"""
	frame_count = os.path.join(script_path, 'orf_frame_count.py')
	output = os.path.join(os.getcwd(), "framed_"+str(os.path.basename(input)))
	if os.path.exists(output) and not overwrite:
		print(f"Prior output found! File: {output}")
		return output
	else:
		frame_cmd = [f"{frame_count}", '-i', f"{input}", '-o', 'framed', '--force']

		subprocess.run(frame_cmd)
	
		assert os.path.exists(output)
		tmp_files.append(output)
		
		return output

def generate_tree(input, overwrite=False):
	"""
	Generate a quick approximate-maximum-likelihood phylogenetic tree for the provided alignment file
	"""
	output = os.path.join(os.getcwd(), str(os.path.basename(input))+".tree")
	if os.path.exists(output) and not overwrite:
		print(f"Prior output found! File: {output}")
		return output
	else:
		tree_cmd = ['fasttree', '-nt', '-gtr', '-gamma', f"{input}"]
		
		with open(output, 'w+') as out:
			subprocess.run(tree_cmd, stdout=out)
		
		assert os.path.exists(output)
		tmp_files.append(output)
		
		return output
	
def hyphy_fit_and_reconstruct(input_seq, input_tree, name, overwrite=False):
	"""
	Use provided HyPhy executable scripts to perform an MG94 fit of the alignemnt with added context from the phylogenetic tree produced by FastTree.
	
	Run provided HyPhy AncestralSequences.bf to perform ancestral reconstruction and parse output JSON which contain a labelled tree and individual node sequences
	"""
	mg94_cmd = None
	hyphy_mg94 = os.path.join(script_path, 'FitMG94.bf')
	hyphy_ancestry = os.path.join(script_path, 'AncestralSequences.bf')
	
	# check if fasta and tree have more than one entry else, exit early
	count_fasta = count_char(f"{input_seq}", ">")
	count_tree = count_char(f"{input_tree}", ":")
	if (count_fasta < 2) and (count_tree < 2):
		return None
	
	# generate MG94 fit
	mg94 = os.path.join(os.getcwd(), f"{name}.mg94.fit")
	if os.path.exists(mg94):
		print(f"Prior output found! File: {mg94}")
		if overwrite:
			os.remove(mg94)
			mg94_cmd = ['hyphy', f"{hyphy_mg94}", '--alignment', f"{input_seq}", '--tree', f"{input_tree}", '--save-fit', f"{mg94}"]
	else:
		mg94_cmd = ['hyphy', f"{hyphy_mg94}", '--alignment', f"{input_seq}", '--tree', f"{input_tree}", '--save-fit', f"{mg94}"]
	
	if mg94_cmd is not None:
		subprocess.run(mg94_cmd)
	if os.path.exists(os.path.join(os.getcwd(), 'errors.log')): # MG94 fit failed (usually due to something wrong with the Newick tree i.e., not enough info)
		os.rename(os.path.join(os.getcwd(), 'errors.log'), os.path.join(os.getcwd(), f"{name}_mg94_errors.log"))
		return None
	
	assert os.path.exists(mg94)
	tmp_files.append(mg94)
	
	# ancestral reconstruction
	ancestor = os.path.join(os.getcwd(), f"{name}.mg94.ancestors.json")
	
	recon_cmd = ['hyphy', f"{hyphy_ancestry}", '--fit', f"{mg94}", '--format', 'JSON', '--output', f"{ancestor}"]
	subprocess.run(recon_cmd)
	if os.path.exists(os.path.join(os.getcwd(), 'errors.log')): # ancestral reconstrution failed (usually due to something wrong with the Newick tree i.e., not enough info)
		os.rename(os.path.join(os.getcwd(), 'errors.log'), os.path.join(os.getcwd(), f"{name}_ancestor_errors.log"))
		with open(f"{ancestor}", 'w+') as blank:
			blank.write(';')

	return ancestor

def clean_tmp(tmp_files):
	for file in tmp_files:
		if os.path.exists(file): os.remove(file)

def run():
	### Parser
	parser = create_parser()
	args = parser.parse_args()
	no_nodes = []
	# correction = os.path.join(script_path, 'orf_correction.py')
	# frame_count = os.path.join(script_path, 'orf_frame_count.py')
	# hyphy_mg94 = os.path.join(script_path, 'FitMG94.bf')
	# hyphy_ancestry = os.path.join(script_path, 'AncestralSequences.bf')
	
	### Parameters
	# check if input is a file or directory
	if os.path.exists(args.input_files) and os.path.isdir(args.input_files):
		fasta_file = [os.path.abspath(file) for file in glob.glob(os.path.join(args.input_files, '*.fa*'))]
		fasta_name = [str(os.path.basename(name)) for name in fasta_file]
	elif os.path.exists(args.input_files) and os.path.isfile(args.input_files):
		fasta_file = [os.path.abspath(args.input_files)]
		fasta_name = [str(os.path.basename(name)) for name in fasta_file]
	else:
		raise IOError("Unknown input")
	
	# check if output exists, change to it to check for intermediates
	if os.path.exists(args.output_dir):
		os.chdir(args.output_dir)
		# if args.overwrite:
			# shutil.rmtree(args.output_dir)
			# os.mkdir(args.output_dir)
			# os.chdir(args.output_dir)
		# else:
			# exit("Output directory exists! Either rename or use '--overwrite'!")
	else:
		os.mkdir(args.output_dir)
		os.chdir(args.output_dir)
	
	for file, name in zip(fasta_file, fasta_name):
		if name.lower().startswith('unclustered'):
			continue
	# remove invalid characters in cluster sequences
	# perform alignment of cluster sequences
		corrected_file = correct_seqs(file, False, args.overwrite)
		align = mafft(corrected_file, name, args.overwrite)
		#align = mafft(file, name, args.overwrite)
		
	# remove invalid characters
		#corrected_mafft = correct_seqs(align, False, args.overwrite)
		
	# trimal automatic trimming
		trimmed_mafft = trim_align(align, args.overwrite)
		#trimmed_mafft = trim_align(corrected_mafft, args.overwrite)
	
	# correct lengths
		length_trimmed = correct_length(trimmed_mafft, args.overwrite)
		
	# remove stop codons
		corrected_trim = correct_seqs(length_trimmed, True, args.overwrite)
		
	# generate approx. ML tree
		fast_tree = generate_tree(trimmed_mafft, args.overwrite)
	
	# Generate MG94 fit and ancestral reconstrction
		ancestor = hyphy_fit_and_reconstruct(corrected_trim, fast_tree, name, args.overwrite)
		if ancestor is None: # mg94 or ancestral reconstruction likely to fail
			no_nodes.append(name)
			continue
		else:
			assert os.path.exists(ancestor)
		
	# parse JSON to extract Newick tree and node sequences
	# json: joint --> ancestral_sequences
	# json: joint --> labeled_tree
		node_seqs = {} # ex. Node1: Node_Seq
		lab_tree = None
		try:
			with open(ancestor) as f:
				nodes = json.load(f)
				lab_tree = nodes['joint']['labeled_tree']
				for n in nodes['joint']['ancestral_sequences']:
					node_seqs[n] = nodes['joint']['ancestral_sequences'][n]
				
			with open(name+".ancestral.tree", 'w+') as outtree:
				outtree.write(str(lab_tree)+';')
				
			with open(name+".ancestral_nodes.fasta", 'w+') as outfasta:
				for n in node_seqs:
					outfasta.write(f">{n}\n{node_seqs[n]}\n")
		except json.decoder.JSONDecodeError: # no output for ancestral reconstruction
			no_nodes.append(name)
			continue
			# with open(name+".ancestral.tree", 'w+') as outtree:
				# outtree.write('')
				
			# with open(name+".ancestral_nodes.fasta", 'w+') as outfasta:
				# outfasta.write('')
			
		
	if len(tmp_files) > 0 and not args.intermed_keep:
		clean_tmp(tmp_files)
		
	os.chdir(start_dir)
	print("DONE")
	
	if len(no_nodes) > 0:
		print(f"No nodes found for:")
		for n in no_nodes:
			print(f"{n}")

if __name__ == '__main__':
	tmp_files = []
	try:
		script_path = os.path.abspath(sys.argv[0].rsplit("/",1)[0])
	except IndexError:
		script_path = os.getcwd() # script in current working directory
	check_scripts(script_path)
	start_dir = os.getcwd()
	run()
	exit()