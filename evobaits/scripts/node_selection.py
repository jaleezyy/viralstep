#!/usr/bin/env python3

### Goal of the script is to take Newick input and run newick_utils
### Using newick_utils: examine distances of tips and nodes from the tree root
### Two general modes: sequenctial distance and median distance
### Sequential distance: Output n FASTA files corresponding to the first closest node to the root, and subsequent FASTA files will include 1 sequence corresponding to the next closest node, repeat until last FASTA containing all node sequences (starting input)
### Median distance: Use maximum tip distance to determine the median (or middle genetic distance), then output node sequences within bottom median (0-median)
### Reverse flag option will pull starting from most distance node in sequential distance mode, and pull the top median (median-max) in median distance mode


import os, sys, shutil
import argparse
import subprocess


def median(values):
	"""
	Take list or dict and determine median
	
	Values will be determined maximum dist using nw_distance -n (looks at distance of each tip to root) + sort
	This will determine the range of genetic distance to work within
	
	"""
	
	if type(values) is dict:
		sorted_values = dict(sorted(values.items(), key=lambda val: val[1]))
		values = [sorted_values[i] for i in sorted_values]
	else:
		values = sorted(values)
	
	mid_a = (len(values) - 1) // 2 # 0-based index
	mid_b = len(values) // 2
	return (values[mid_a] + values[mid_b]) / 2
	
def dist_to_dict(std):
	"""
	Will use sorted output from nw_distance -s i -n which will provide a list of nodes and their corresponding distance from the root
	
	Expected output will contain: label\tdistance
	"""
	dist_dict = {}
	for line in std.stdout.decode("utf-8").split('\n'):
		parse = line.split("\t")
		if len(parse) == 2:
			label = parse[0].strip()
			dist = parse[1].strip()
			if label not in dist_dict:
				dist_dict[label] = float(dist) # Newick int only?
		else:
			continue
	
	return dist_dict

def out_indiv(seqs, node_dict, outdir, rev=False):
	node_dict = dict(sorted(node_dict.items(), key=lambda val: val[1], reverse=rev))
	node_list = list(node_dict.keys())
	num_seq = len([1 for line in open(seqs) if line.startswith(">") and "root" not in line.lower().strip(">")])
	if num_seq < 1:
		print("No sequences found!")
		sys.exit(1)
	seq_collection = {}
	
	with open(seqs) as infasta:
		while True:
			line1 = infasta.readline() #>header
			line2 = infasta.readline()
			if not line2: break # EOF
			assert line1.startswith(">")
			if line1.strip().strip(">").lower() == 'root': continue
			seq_collection[line1.strip().strip(">")] = line2.strip()
	#try:
	assert len(seq_collection) == num_seq # each header should be unique
	#except AssertionError:
	#	print(seq_collection)
	#	print(f"Number of sequences parsed: {len(seq_collection)}")
	#	print(f"Number of sequences in file: {num_seq}")
		
	for i in range(1,num_seq+1):
		output = os.path.join(outdir, f"{i}_distance.fasta")
		subset = node_list[0:i] # 0-based indexing
		with open(output, 'w+') as outfasta:
			for header in subset:
				outfasta.write(f">{header}\n{seq_collection[header]}\n")
				
	print(node_dict)
	
def out_median(seqs, node_dict, median, outdir, rev=False):
	if median == 0:
		print("No median distance generated! Cannot produce output!")
		sys.exit(1)
	
	num_seq = len([1 for line in open(seqs) if line.startswith(">") and "root" not in line.lower().strip(">")])
	if num_seq < 1:
		print("No sequences found!")
		sys.exit(1)
	
	### TODO: optimize by using sorted list to minimize parsing the whole FASTA file
	node_dict = dict(sorted(node_dict.items(), key=lambda val: val[1], reverse=rev))
	output = os.path.join(outdir, 'median_distance.fasta')
	
	with open(seqs) as infasta, open(output, 'w+') as outfasta:
		while True:
			line1 = infasta.readline() #>header\n
			line2 = infasta.readline() #seq\n
			if not line2: break # EOF
			assert line1.startswith(">")
			if line1.strip().strip(">").lower() == 'root': continue
			
			if rev:
				if float(node_dict[line1.strip().strip(">")]) >= median:
					outfasta.write(f"{line1}{line2}")
			else:
				if float(node_dict[line1.strip().strip(">")]) <= median:
					outfasta.write(f"{line1}{line2}") 
	
	print(node_dict)
	print(f"Median: {median}")
	
def check_nw():
	params = ['nw_distance', '-h']
	try:
		subprocess.run(params, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError:
		print("Newick_utils found however, something went wrong! Check installation!")
		sys.exit(1)
	except FileNotFoundError:
		print("Missing Newick_utils installation!")
		sys.exit(1)
	pass

def run(args):
	### Parameters
	fasta_name = os.path.basename(args.input_fasta_file)
	tip_params = ['nw_distance', '-n', str(args.input_tree)] # tip distances
	node_params = ['nw_distance', '-n', '-s', 'i', str(args.input_tree)] # node distances
	
	### Check for Newick_utils *critical*
	check_nw()
	
	### 
	try:
		os.mkdir(args.output_dir)
	except FileExistsError:
		if args.overwrite:
			shutil.rmtree(args.output_dir)
			os.mkdir(args.output_dir)
		else:
			print(f"Output directory '{args.output_dir}' exists! Use --overwrite or rename, and try again!")
			sys.exit(1)
	
	### Pull node distances (and sort)
	node = subprocess.run(node_params, stdout=subprocess.PIPE)
	node_distances = dist_to_dict(node)
	
	### Check mode (if median, calculate)
	if args.mode == "median":
		tip = subprocess.run(tip_params, stdout=subprocess.PIPE)
		tip_distances = dist_to_dict(tip)
		median_dist = median(tip_distances)
		print(tip_distances)
	elif args.mode == 'node_median':
		median_dist = median(node_distances)
	else: # default, seq
		median_dist = 0
	
	### Output based on mode
	if args.mode == "median" or args.mode == "node_median":
		out_median(args.input_fasta_file, node_distances, median_dist, args.output_dir, args.rev)
	else: # default, seq
		out_indiv(args.input_fasta_file, node_distances, args.output_dir, args.rev)

def create_parser():
	parser = argparse.ArgumentParser(description='Run Newick_utils on a given tree alongside node sequences to parse and extract seqeunces for bait design.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file containing node sequences")
	parser.add_argument('-t', '--tree', dest='input_tree', required=True, help='Input Newick tree file')
	parser.add_argument('-o', '--outdir', dest='output_dir', default='node_parse', help="Output directory. Default='node_parse'")
	parser.add_argument('--mode', choices=['seq', 'median', 'node_median'], default='seq', help="Parsing mode. 'Seq' = Sequenctial mode. Output FASTA files will contain a single sequence that is closest to the tree root. Additional FASTA output will include one additional sequence until end of tree. If the '--rev' flag is used, output will begin from the furthest node from the root. 'Median' = Median mode. The median distance from the root will be determined from the tips of the tree. Output FASTA will contain all node sequences closest (or furthest if '--rev') to the median. Node Median = Same principle as Median mode, but median will instead be calculated using the distances of the nodes rather than tips")
	parser.add_argument('--rev', action='store_true', help="Reverse order of sequences pulled")
	parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output directory!')
	
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
	run(args)