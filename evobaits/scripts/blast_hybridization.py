#!/bin/env python3

### OUTPUT SHOULD BE COMPATIBLE WITH "viralstep plots"
### TODO: Improve logging

from Bio import SeqIO
import argparse
import os, sys, subprocess
import glob, shutil
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot

def create_parser():
	parser = argparse.ArgumentParser(prog=".", description='Calculate coverage using BLAST local alignment to reflect potential hybridization outcomes')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file")
	parser.add_argument('-r', '--reference', dest='reference_file', required=True, help="Reference sequences for hybridization targets")
	parser.add_argument('-o', '--output', dest="output_file", default='merged_coverage.txt', required=True, help="Output file(s)")
	parser.add_argument('-d', '--database', dest='reference_db', required=False, help="ViralSTEP DB file containing taxonomy information for reference sequences")
	parser.add_argument('-p', '--proportion', dest='cov_proportion', default=70.0, help='Percent of alignment between bait and reference sequence to consider it as successful hybridization. Default = 70.0')
	parser.add_argument('-t', '--threads', dest="num_threads", required=False, default=1, type=int, help="Number of threads for BLASTn. Default = 1")
	parser.add_argument('-w', '--word-size', dest='word_size', default=11, help='BLASTn parameter. Word size used to align sequences. Default BLASTn = 11') 
	parser.add_argument('--alt-stats', dest='additional_stats', action='store_true', help='Output additional stats files, replacing ViralSTEP-specific stats')
	parser.add_argument('--keep-intermediate', dest="output_pre", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	parser.add_argument('--overwrite', action='store_true', help='Overwrite any existing intermediate files and databases. By default, if these files are found, the file generation steps are otherwise skipped')
	parser.add_argument('--keep-zeros', dest='keep_zeros', action='store_true', help="Add zero coverage statistics into final output.")
	parser.add_argument('--skip-plots', dest='skip_plots', action='store_true', help="FOR DEBUGGING PURPOSES: disable per sequence coverage plot generation.")
	return parser

def count_seq(input):
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq

def median(values):
	"""
	Take list or dict and determine median
	
	"""
	
	if type(values) is dict:
		sorted_values = dict(sorted(values.items(), key=lambda val: val[1]))
		values = [sorted_values[i] for i in sorted_values]
	else:
		values = sorted(values)
	
	mid_a = (len(values) - 1) // 2 # 0-based index
	mid_b = len(values) // 2
	return (values[mid_a] + values[mid_b]) / 2

def load_reference_stats(file):
	"""
	Store accession: seq_length for each sequence for a given file
	"""
	
	ref = {}
	
	for seq_record in SeqIO.parse(file, 'fasta'):
		id = seq_record.id
		length = len(seq_record.seq)
		
		ref[id] = length
		
	return ref

def load_reference(file, subset=None):
	"""
	Output defaultdict object where accession ID stores each nucleotide position + coverage.
	
	If a list of accessions is provided, defaultdict object will only contain those keys
	"""
	
	ref = defaultdict(dict) # {acc: {pos: cov}}
	
	for seq_record in SeqIO.parse(file, 'fasta'):
		id = seq_record.id
		length = len(seq_record.seq) # length of sequence
		if subset is not None and isinstance(subset, list):
			if id in subset:
				for pos in range(1,length+1):
					ref[id][pos] = 0 # set coverage to 0
			else:
				continue
		else:
			for pos in range(1,length+1):
				ref[id][pos] = 0 # set coverage to 0
			
	return ref

def load_taxonomy(refdb):
	"""
	Given a ViralSTEP DB file, pull accession ID and corresponding virus family.
	
	Breakdown of ViralSTEP DB: AccID, TaxID, Baltimore, Virus Family. Minimum is AccID, Virus Family.
	"""
	virus_families = {}
	with open(refdb) as db:
		for line in db:
			row = line.split(',')
			acc = row[0]
			# tax = row[1]
			# balt = row[2]
			fam = row[-1] # row[3]
			virus_families[acc] = fam
			
	return virus_families
			

def makeselfblastdb(input, ref_name, output_loc):
	start_dir = os.getcwd()
	fasta_file = os.path.basename(input)
	db_path = os.path.join(output_loc, f"{ref_name}_blastdb") # out_dir is an absolute path
	if os.path.isdir(db_path):
		print(f"Updating existing {ref_name} DB")
		shutil.rmtree(db_path) # delete previous
		os.mkdir(db_path)
		os.chdir(db_path)
	else:
		os.mkdir(db_path)
		os.chdir(db_path)
	print(f"Generating BLASTDB of provided reference FASTA file: {fasta_file}")
	dbline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=input, out=f"{ref_name}")
	stdout, stderr = dbline()
	os.chdir(start_dir)
	return db_path

def blastn(input, db, ref_name, output_loc, threads, wsize=11, ref_taxonomy=None):
	fasta_name = os.path.basename(input)
	abs_out = os.path.abspath(output_loc) # absolute path likely not necessary
	blast = os.path.join(abs_out, f'{ref_name}.all.blast.out')
	
	print("Running BLASTn")
	### std: 'qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore'
	if ref_taxonomy is not None:
		fmt = '10 std staxids sscinames sblastnames sskingdoms'
		index = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames', 'sblastnames', 'sskingdoms']
	else:
		fmt = '10 std'
		index = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
	#blastline = NcbiblastnCommandline(query=input, db=f"{os.path.join(db, ref_name)}", outfmt=f'{fmt}', out=blast, num_threads=threads, evalue=0.05, max_target_seqs=100000, word_size=wsize)
	blastline = NcbiblastnCommandline(query=input, db=f"{os.path.join(db, ref_name)}", outfmt=f'{fmt}', out=blast, num_threads=threads, max_target_seqs=100000, task='blastn-short', evalue=0.01)
	stdout, stderr = blastline()
	
	return blast, index

def coverage(ref, ref_length, bait_length, hit, cov=70.0):
	"""
	Using the HSP data (hit), we can infer the proportion (using ref_length and bait_length) of the bait that hybridizes/aligns (query; qseqid) and the nucleotide positions covered on the reference (subject; sseqid), returning the updated coverage defaultdict (ref). 
	"""
	# check evalue to ensure quality hit
	if float(hit['evalue']) >= 0.01:
		return ref
		
	# determine proportion
	if bait_length > 0:
		prop = (hit['pident']*hit['length'])/bait_length
	else:
		prop = hit['pident']
	
	# if proportion sufficient, find boundaries on subject
	# we'll assume sstart and send have full coverage regardless of similairy + align length
	if prop < float(cov):
		return ref
	else:
		# update coverage
		for pos in range(hit['sstart'], hit['send']+1):
			ref[hit['sseqid']][pos]+=1
	
	# if sskingdoms found, add taxonomy information (TODO)
	
	return ref

def coverage_calc(ref_cov, ref_lengths, keep_zeros=False):
	"""
	Input is a defaultdict containing {acc: {pos: cov}}. 'ref_length' provides all reference accessions.
	
	Per accession, store a tuple with (depth, breadth) stored in a collective dict object.
	
	All other accessions without hits will be assumed 0 (fold or percent) for depth and breadth of coverage, respectively.
	"""
	dep_bread_cov = {} # {acc: (depth, breadth, etc...)}
	cov = []
	for acc in ref_lengths:
		if acc in ref_cov:
			count = 0 # all non-zero positions
			depth_count = 0
			genome_len = ref_lengths[acc]
			for pos in ref_cov[acc]:
				value = ref_cov[acc][pos] # value at position
				if value > 0:
					count += 1
					depth_count = depth_count + value
					cov.append(value)
			breadth = float((count/genome_len)*100)
			depth = float(depth_count/genome_len)
			try:
				depth_covered = float(depth_count/count)
			except ZeroDivisionError:
				depth_covered = 0
			
			med = median(cov)
			median_covered = median([c for c in cov if c != 0]) # median for non-zero coverage
			
			dep_bread_cov[acc] = (depth, breadth, depth_covered, med, median_covered)
		else: # no hit; a zero
			if keep_zeros:
				dep_bread_cov[acc] = (0,0,0,0,0)
			else:
				continue
	return dep_bread_cov
	

def coverage_plot(dbcov, output, name, keep_zeros=False, skip=False, overwrite=False):
	"""
	Input is a defaultdict containing {acc: {pos: cov}}. 'ref_length' provides all reference accessions.
	
	Per accession with a hit, generate a matplotlib plot. Optionally can generate plots for 0-fold depth or 0-percent breadth of coverage.
	"""
	out_dir = os.path.join(output, f"{name}_plots")
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	elif overwrite:
		shutil.rmtree(out_dir)
		os.mkdir(out_dir)
	else: # if found, but no overwrite flag: skip
		print("WARNING: Plots directory already exists! Use '--overwrite' to replace!")
		return 0
	
	if not skip:
		for acc in dbcov:
			cleaned_acc = '_'.join(str(acc).split('_')[:2])
			pos = [i for i in dbcov[acc]]
			depth = [dbcov[acc][j] for j in dbcov[acc]]
			assert len(pos) == len(depth)
			
			print(f"Generating coverage plot for {cleaned_acc}!")
			cov_plot = matplotlib.pyplot.figure()
			matplotlib.pyplot.plot(pos, depth)
			cov_plot.suptitle(f"{cleaned_acc}", fontsize=16)
			matplotlib.pyplot.xlabel("Nucleotide Position", fontsize=12)
			matplotlib.pyplot.ylabel("Depth of Coverage", fontsize=12)
			
			cov_plot.savefig(os.path.join(out_dir, f"{cleaned_acc}_coverage_plot.png"), dpi=300)
	else:
		print("Skipping depth of coverage plot!")
		
	
def coverage_output(dbcov, output, ref_taxonomy=None, overwrite=False, alt_stats=False):
	"""
	Final output should be found where 'plots' directory is. Parameter 'output' refers to args.output_file which will have path & file. 
	
	Breakdown (tab-delimited):
		Accession, Genome Completeness, Average Depth of Coverage, Bait Origin (default = ancestral), Hit Origin(s).
		
	Bait Origin refers to virus family of reference sequence. Legacy naming from ViralSTEP; may rename.
	
	Alternate stats breakdown (tab-delimited):
		Accession, Genome Completeness, Average Depth of Coverage of Covered Regions, Bait Origin (default = ancestral), Median, Median of Covered Regions
	"""
	add_out = f"additional_stats_for_covered_regions_{output}"
	if os.path.exists(output) and not overwrite:
		print(f"WARNING: Output file {output} already exists! Use '--overwrite' to replace!")
		return 0
	elif os.path.exists(output) and overwrite:
		os.remove(output)

	if alt_stats:
		if os.path.exists(add_out) and not overwrite:
			print(f"WARNING: Output file {add_out} already exists! Use '--overwrite' to replace!")
			return 0
		elif os.path.exists(add_out) and overwrite:
			os.remove(add_out)
		out2 = open(add_out, 'w+')
		out2.write(f"Accession\tGenome Completeness\tAverage Depth of Coverage against Covered Regions\tBait Origin\tMedian\tMedian of Covered Regions\n")
	
	with open(output, 'w+') as out:
		print("Generating summary file of coverages!")
		out.write(f"Accession\tGenome Completeness\tAverage Depth of Coverage\tBait Origin\tHit Origin(s)\n")
		for acc in dbcov:
			cleaned_acc = '_'.join(str(acc).split('_')[:2])
			depth = format(dbcov[acc][0], '.5f')
			breadth = format(dbcov[acc][1], '.5f')
			
			if ref_taxonomy is not None:
				tax = str(ref_taxonomy[acc]).strip()
			else:
				tax = 'Ancestral'
			
			out.write(f"{cleaned_acc}\t{breadth}%\t{depth}x\t{tax}\t{{}}\n")
			
			if alt_stats:
				depth_covered = format(dbcov[acc][2], '.5f')
				med = format(dbcov[acc][3], '.5f')
				median_covered = format(dbcov[acc][4], '.5f')
				
				out2.write(f"{cleaned_acc}\t{breadth}%\t{depth_covered}x\t{tax}\t{med}x\t{median_covered}x\n")
	
	if alt_stats:
		out2.close()

def clean_tmp(tmp, delete=True):
	pass

def run(args):
	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	ref_file = os.path.abspath(args.reference_file)
	ref_name = str(os.path.basename(args.reference_file))
	num_seq = count_seq(fasta_file)
	try:
		out_dir = os.path.abspath(args.output_file).rsplit("/", 1)[0]
	except IndexError:
		out_dir = os.getcwd()
	threads = int(args.num_threads)
	
	### Parse reference sequence to collect seqids
	### acc: seq_length
	ref_stats = load_reference_stats(ref_file)
	
	
	### Count bait length, should be consistent, will change mode of hybridization calculation
	### Mode 'prop' calculates proportion of bait that aligns (similar to bait redundancy) if bait_length > 0
	### Mode 'sim' uses alignment sequence similarity to estimate
	bait_seq_sizes = list(set([len(seq_record.seq) for seq_record in SeqIO.parse(fasta_file, 'fasta')]))
	if len(bait_seq_sizes) == 1:
		bait_length = bait_seq_sizes[0]
		# prop will look at [prop]ortion of bait that aligns within local alignment to judge hybridization
	else:
		bait_length = 0
		# sim will use only sequence [sim]ilarity for each BLAST hit i.e., local alignment, to gauge hybridization
	
	### If ViralSTEP DB file found, generate separate parameter 
	### Store virus family information
	if args.reference_db is not None:
		ref_taxonomy = load_taxonomy(args.reference_db)
	else:
		ref_taxonomy = None
	#virus_families = {} ### TODO
	
	### Make BLASTDB of reference sequences
	### BLASTDB found in ref_name/ref_name.n*
	expected_db = [f'{ref_name}.ndb', f'{ref_name}.nhr', f'{ref_name}.nin', f'{ref_name}.not', f'{ref_name}.nsq', f'{ref_name}.ntf', f'{ref_name}.nto']
	expected_blast = os.path.join(out_dir, f'{ref_name}.all.blast.out')
	if all([os.path.exists(os.path.join(ref_name, i)) for i in expected_db]) and not args.overwrite:
		# if found, do not regenerate, skip step
		print("BLASTDB found and will not be overwritten!")
	else:
		# run makeblastdb (allow overwrite step)
		db_path = makeselfblastdb(ref_file, ref_name, out_dir)
		assert all([os.path.exists(os.path.join(db_path, f)) for f in expected_db])
	
	if os.path.exists(expected_blast) and not args.overwrite:
		# if found, do not regenerate, skip step
		print("BLAST outpui found! Will parse existing file!")
		blast_out = expected_blast
		# check first line for commas: 11 = 12 columns; 15 = 16 columns
		with open(blast_out) as f:
			first_line = f.readline()
			if first_line.count(',') == 11:
				index = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
			elif first_line.count(',') == 15:
				index = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames', 'sblastnames', 'sskingdoms']
			else:
				index = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
	else:
		# run BLASTN
		if args.overwrite and os.path.exists(expected_blast):
			print("Overwriting existing BLAST results!")
			os.remove(expected_blast)
		blast_out, index = blastn(fasta_file, db_path, ref_name, out_dir, threads, args.word_size, ref_taxonomy)
		assert os.path.exists(blast_out)
		assert f"{os.path.basename(blast_out)}" == f"{os.path.basename(expected_blast)}"
		
	### Parse BLAST output (comma-separated)
	blast_df = pd.read_csv(blast_out, header=0, names=index)
	print(blast_df)
	
	hits = list(blast_df.drop_duplicates(subset=['sseqid'], keep='first')['sseqid'])
	#print(hits)
	ref_cov = load_reference(ref_file, hits)
	#print(ref_cov)
	for h in hits:
		#print(h)
		subset = blast_df[blast_df['sseqid'] == h]
		#print(subset)
		subset.apply(lambda row: coverage(ref_cov, ref_stats, bait_length, row, args.cov_proportion), axis=1)
	
	# with open('test.log', 'w+') as out:
		# print(ref_cov, file=out)
		
	### Use ref_cov to determine stats
	### Breadth/Depth of coverage
	depth_breadth = coverage_calc(ref_cov, ref_stats, args.keep_zeros) # acc: (depth, bredth)
	
	### Use ref_cov to generate plots per reference (i.e., subject)
	coverage_plot(ref_cov, out_dir, fasta_name, args.keep_zeros, args.skip_plots, args.overwrite)
	
	### Final output file containining coverage values per accession
	coverage_output(depth_breadth, args.output_file, ref_taxonomy, args.overwrite, args.additional_stats)

	

if __name__ == '__main__':
	try:
		script_path = sys.argv[0].rsplit("/",1)[0]
	except IndexError:
		script_path = os.getcwd() # script in current working directory
	start_dir = os.getcwd()
	tmp_files = []
	### Parser
	parser = create_parser()
	args = parser.parse_args()
	run(args)
	sys.exit(0)