#!/bin/env python

import argparse
import os, sys, shutil, glob
import datetime, time
from Bio import SeqIO
import urllib.request
from collections import defaultdict
import subprocess
import pandas
import pysam
import matplotlib.pyplot
import pybedtools

# custom modules
try:
	from scripts.fasta2bed import fasta2bed 
	from scripts.split_fasta import split_fasta, batch_iterator
except ModuleNotFoundError: # running viralstep.py directly
	from fasta2bed import fasta2bed 
	from split_fasta import split_fasta, batch_iterator


### Goal is to take a bait set and examine coverage across multple virus species
### Align targets will be determined using the references comprising the bait set
### Will download reference sequences to one or more files and generate alignment db (HiSAT2 or Bowtie2)
### samtools mpileup (per-base) 
### Compute genome completeness (% of bases with >=1 bait covered)
### Compute genome depth of coverage (average, median) (# of baits per base)

### TODO: Add a version check for software to ensure it exists prior to running
### TODO: Add checks for aligner (i.e., Bowtie2)


def create_parser():
	parser = argparse.ArgumentParser(prog='viralstep coverage', description='Calculate bait coverage statistics against reference sequences representatitve of the given bait set.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file")
	#parser.add_argument('-p', '--proportion', dest="input_prop_file", required=False, help="Input file containing breakdown of bait set by virus family. If not provided, classification pipeline will be run using database file")
	parser.add_argument('-d', '--database', dest="input_db_file", required=True, help="Input CSV file containing individual bait accessionID, taxID, baltimore classification, virus family. Must align with reference FASTA if --reference flag provided. Otherwise must align with bait FASTA.")
	parser.add_argument('-r', '--reference', dest="manual_ref", required=False, default=None, help="Single FASTA file containing all desired reference sequences. If not provided, then bait headers will determine reference sequences.")
	parser.add_argument('-a', '--aligner', dest="align_tool", required=False, default="bowtie2", help="Specify aligner to use between Bowtie2/HiSAT2/BWA/KMA. Default=bowtie2")
	#parser.add_argument('-o', '--out', dest="output_prefix", required=True, help="Output prefix for directories. Default=FASTA filename.")
	parser.add_argument('-t', '--threads', dest="num_threads", required=False, default=1, help="Number of threads for alignment. Default = 1")
	parser.add_argument('--ngscat', dest="run_ngscat", action="store_true", required=False, help="Run ngsCAT analysis. Slows down overall runtime.")
	parser.add_argument('--coveragethrs', dest="coveragethresholds", required=False, default='1,5,10,20,30', help="Optional argument for ngsCAT. Comma separated list of real numbers (do not leave spaces between) indicating coverage thresholds to be used when calculating percentages of covered bases (first graph in the ngsCAT report). Default=1,5,10,20,30.")
	parser.add_argument('--keep-intermediate', dest="output_pre", action="store_true", required=False, help="Keep intermediate temporary files") ### will delete at end of program otherwise
	parser.add_argument('-k', '--k-mer', dest="k_size", required=False, default="16", help="Optional argument for KMA aligner. Will index and align accordingly. Default k-mer size of 16 will be used if not specified.")
	parser.add_argument('--skip-plots', dest='skip_plots', action='store_true', help="FOR DEBUGGING PURPOSES: disable per sequence coverage plot generation.")
	return parser

def count_seq(input):
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq
	
def download_ref(accession, ref_dir):
	os.chdir(ref_dir)
	# individual references
	for acc in accession:
		with urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=txt" %(acc)) as response, open("%s.fasta" %(acc), 'wb') as out_file:
			shutil.copyfileobj(response, out_file)
			print(acc)
			time.sleep(0.3)
	# combined references
	print("\nCreating merged reference FASTA\n")
	with open("merged.fasta", 'w') as merged:
		for file in glob.glob("*.fasta"):
			if "merged.fasta" == file: continue
			for line in open(file):
				merged.write(line)
	
	os.chdir(start_dir)

def split_ref(ref_fasta, ref_dir):
	os.chdir(ref_dir)
	# copy over ref_fasta
	shutil.copy(ref_fasta, os.path.join(os.getcwd(), "merged.fasta"))
	print("Splitting given reference FASTA\n")
	split_fasta("merged.fasta", 1) # should split into group0#_merged.fasta
	for file in glob.glob(os.path.join(os.getcwd(), "group_0*_merged.fasta")):
		if "merged.fasta" == file: continue
		num_seq = count_seq(file)
		if num_seq == 1:
			for seq_record in SeqIO.parse(file, 'fasta'):
				if str(seq_record.id).count("_") > 1:
					name = str(seq_record.id).rsplit("_",1)[0] # assumes bait headers Acc_ID_1-80
				else: 
					name = str(seq_record.id)
		try:
			os.replace(file, '%s.fasta' %(name))
		except OSError: # file with accession already exists; if no exception, then silently replaced
			os.remove(file) 
	
	os.chdir(start_dir)

def index_ref(ref_dir, k_size, index):
	os.chdir(ref_dir)
	for file in glob.glob("*.fasta"):
		name = str(os.path.basename(file)).rsplit(".",1)[0]
		if "hisat2" == index.lower():
			exit("Add support for HiSAT2")
			#subprocess.run(['hisat2-build', '%s' %(file), '%s' %(str(file).rsplit(".",1)[0])])
		elif "bwa" == index.lower():
			exit("Add support for BWA")
			#subprocess.run(['bwa', 'index', '-p', '%s' %(str(file).rsplit(".",1)[0]), '%s' %(file)])
		elif "kma" == index.lower():
			subprocess.run(['kma', 'index', '-i', '%s' %(file), '-o', '%s' %(name), '-k', '%s' %(k_size)])
		elif "bowtie2" == index.lower():
			subprocess.run(['bowtie2-build', '%s' %(file), '%s' %(name)])
		else:
			print("Unknown alignment tool. Running Bowtie2-build")
			subprocess.run(['bowtie2-build', '%s' %(file), '%s' %(name)])
	os.chdir(start_dir)

# Analyze BWA, Bowtie2, HiSAT2 alignment
def completeness(bam, acc, cov_dir, bait_hits, bait_origin, skip=False):
	bam_file = pybedtools.BedTool(bam)
	genomecov = bam_file.genome_coverage(d=True)
	### TODO: Add ability to generate plot using 'pos' and 'depth'
	df = genomecov.to_dataframe(names=['chrom', 'pos','depth'])
	
	if not skip:
		print("\nGenerating depth of coverage plot for %s" %(acc))
		cov_plot = matplotlib.pyplot.figure()
		matplotlib.pyplot.plot(df['pos'], df['depth'])
		cov_plot.suptitle("%s" %(acc), fontsize=16)
		matplotlib.pyplot.xlabel("Nucleotide Position", fontsize=12)
		matplotlib.pyplot.ylabel("Depth of Coverage", fontsize=12)
		cov_plot.savefig(os.path.join(cov_dir, "plots", "%s_coverage_plot.png" %(acc)), dpi=300)
	else:
		print("Skipping depth of coverage plot!")
	
	genome_len = len(df)
	count = 0
	depth_count = 0
	
	# determine bait cross-hybridization
	# bait_origin = acc --> virus family
	print("Determining average coverages and taxonomic identification between bait(s) and target\n")
	bait_hit_family = str(bait_hits[acc])
	
	with open(os.path.join(cov_dir, '%s.txt' %(acc)), 'w+') as out:
		out.write("Accession\tGenome Completeness\tAverage Depth of Coverage\tAverage Depth of Coverage Across Non-Zero\tBait Origin\tHit Origin(s)\n")
		for value in df['depth']:
			if value > 0:
				count += 1
				depth_count = depth_count + value
		complete = float((count/genome_len)*100)
		avg_depth = float(depth_count/genome_len)
		try:
			avg_depth_non_zero = float(depth_count/count)
		except ZeroDivisionError:
			avg_depth_non_zero = 0.00
		# count represents the length of the genome where depth > 0
		# avg_depth_non_zero represent DoC across non-zero regions (i.e., to account for gapped coverage)
		
		out.write("%s\t%s%%\t%sx\t%sx\t%s\t%s" %(acc, format(complete, '.2f'), format(avg_depth, '.2f'), format(avg_depth_non_zero, '.2f'), bait_origin, bait_hit_family))
		


def align_baits(baits, ref_dir, cov_dir, tool, virfam_db, threads, skip=False):
	warn = False
	# generate virus family acc --> virus family
	accID2family = {} # ex. {MN908947.3: coronaviridae}
	for line in open(virfam_db):
		col = line.split(",")
		if col[0].count("_") > 1:
			accID = col[0].rsplit("_",1)[0] # assume bait headers ex. Acc_ID_1-80
		else:
			accID = col[0]
		accID2family[accID] = col[-1].strip("\n")
	
	os.chdir(ref_dir)
	for file in glob.glob("*.fasta"):
		bait_hits = defaultdict(dict) # {virus_family_ref = {virus_family_baits: count}}
		#db_path = str(file).split(".",1)[0]
		db = str(os.path.basename(file)).rsplit(".",1)[0] # reference
		if "hisat2" == tool.lower():
			exit("Add support for HiSAT2")
			#subprocess.run(['hisat2', '%s' %(file), '%s' %(str(file).rsplit(".",1)[0])])
		elif "bwa" == tool.lower():
			exit("Add support for BWA")
			#subprocess.run(['bwa', 'mem', '-p', '%s' %(str(file).rsplit(".",1)[0]), '%s' %(file)])
		elif "kma" == tool.lower():
			with open('%s.sam' %(db), 'w+') as kma_out:
				subprocess.run(['kma', '-sam', '-nc', '-nf', '-mem_mode', \
									'-o', '%s' %(db), '-t_db', '%s' %(db), '-i', '%s' %(baits), \
									'-t', '%s' %(threads), '-a'], stdout = kma_out)
		elif "bowtie2" == tool.lower():
			subprocess.run(['bowtie2','-x','%s' %(db), '-f', '%s' %(baits), '-S', '%s.sam' %(db), \
								'-p', '%s' %(threads), '--end-to-end', '-N', '1', '-L 32', '-a'])
		else:
			print("Unknown alignment tool. Running Bowtie2")
			subprocess.run(['bowtie2','-x','%s' %(db), '-f', '%s' %(baits), '-S', '%s.sam' %(db), \
								'-p', '%s' %(threads), '--end-to-end', '-N', '1', '-L 32', '-a'])
			#subprocess.run(['bowtie2','-x','%s' %(db), '-f', '%s' %(baits), '-S', '%s.sam' %(db), \
								#-p', '%s' %(threads), '--local', '-N', '1', '-L 32', '-a'])
		
		pysam.samtools.view('-b', '-o', '%s.bam' %(db), '-S', '%s.sam' %(db), catch_stdout=False)
		pysam.sort('-o', '%s.bam' %(db), '%s.bam' %(db))
		pysam.index('%s.bam' %(db))
		os.remove('%s.sam' %(db)) # remove unneeded SAM file
		
		# bait hybridization
		if db != "merged":
			try:
				ref_origin = accID2family[db]
			except KeyError:
				try:
					ref_origin = accID2family[db.split(".")[0]] # to account for accession.version (from manual reference input)
				except KeyError:
					print("\nMissing accession %s within provided database file! Skipping analysis!" %(db))
					warn = True
					continue
			for line in pysam.samtools.view('%s.bam' %(db), '-F', '4', catch_stdout=True).split("\n")[:-1]:
				if line.split("\t")[0].count("_") > 1:
					hit = line.split("\t")[0].rsplit("_",1)[0]
				else:
					hit = line.split("\t")[0]
				try:
					if accID2family[hit] in bait_hits[db]:
						bait_hits[db][accID2family[hit]] = bait_hits[db][accID2family[hit]] + 1
					else:
						bait_hits[db][accID2family[hit]] = 1
				except KeyError: # rethink logic here
					if "Unknown" in bait_hits[db]:
						bait_hits[db]["Unknown"] = bait_hits[db]["Unknown"] + 1
					else:
						bait_hits[db]["Unknown"] = 1
					#print("Missing hit accession %s in database file! Skipping count!" %(hit))
					warn = True
			
		# bedtools coverage
		if db != "merged": completeness('%s.bam' %(db), db, cov_dir, bait_hits, ref_origin, skip)
		
	os.chdir(start_dir)
	return warn

### TODO: Add breakdown statistics
def prop_calc(input_fasta, input_db): # proportion of given bait set
	prop_missing = True
	pass
	return prop_missing

# def stats_calc(baits, pileup): # average / median / ngscat
	# pass
	
def ngscat(ref_dir, cov_dir, cov_threshold, threads):
	# prep bed files
	os.chdir(ref_dir)
	for file in glob.glob('*.fasta'):
		fasta2bed(file)
	os.chdir(cov_dir)
	for file in glob.glob(os.path.join(ref_dir, '*.bam')):
		acc = str(os.path.basename(file)).rsplit(".", 1)[0]
		subprocess.run(['python', '%s' %(os.path.join(directory, 'ngscat', 'ngscat.py')),\
						'--bams', '%s' %(os.path.abspath(file)), \
						'--bed', '%s' %(os.path.abspath(os.path.join(ref_dir, acc+'.bed'))), \
						'--out', '%s' %(os.path.abspath(os.path.join(cov_dir, acc))), \
						'--reference', '%s' %(os.path.abspath(os.path.join(ref_dir, acc+'.fasta'))), \
						'--coveragethrs', '%s' %(cov_threshold), \
						'--threads', '%s' %(threads) \
						])
					

def clean(fasta_name, ref_dir):
	for file in glob.glob(os.path.join(ref_dir, '*.bam')):
		filename = str(os.path.basename(file))
		if os.path.exists(os.path.join('/tmp', filename)): 
			os.remove(os.path.join('/tmp', filename))
	if os.path.exists(ref_dir): 
		shutil.rmtree(ref_dir)

def run():
	### Parser and parameters
	prop_missing = False
	parser = create_parser()
	args = parser.parse_args()
	
	# ngscat requirements
	if args.run_ngscat is True:
		try:
			import sets
			import xlwt
			import progressbar
		except ModuleNotFoundError:
			#print("Missing modules for ngsCAT, auto-install?")
			auto_install = input("Missing modules for ngsCAT, auto-install? [y/n]: ")
			if str(auto_install).lower() in ("y", "yes"):
				subprocess.run(['python', '-m', 'pip', 'install', 'pysam', 'sets', 'xlwt', 'matplotlib', 'progressbar', 'pybedtools'])
			else:
				exit("""
		Please ensure the following modules are installed:
			pysam
			sets
			xlwt
			matplotlib
			progressbar
			pybedtools
			""")
	
	threads = int(args.num_threads)
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = str(os.path.basename(fasta_file))
	
	ref_dir = os.path.join(start_dir, fasta_name+"_coverageref")
	cov_dir = os.path.join(start_dir, fasta_name+"_coveragealign")
	
	# directories 
	if os.path.exists(ref_dir):
		shutil.rmtree(ref_dir)
		os.mkdir(ref_dir)
	else: 
		os.mkdir(ref_dir)
	if os.path.exists(cov_dir):
		shutil.rmtree(cov_dir)
		os.makedirs(os.path.join(cov_dir, "plots"))
	else: 
		os.makedirs(os.path.join(cov_dir, "plots"))
	
	# reference genome sequences (either by file or by download)
	if args.manual_ref is not None:
		split_ref(os.path.abspath(args.manual_ref), ref_dir)
	else:
		ref_accessions = set()
		for seq_record in SeqIO.parse(fasta_file, 'fasta'):
			ref_accessions.add(seq_record.id.rsplit("_",1)[0])
		# Download reference FASTAs to map
		print("Downloading references\n")
		download_ref(ref_accessions, ref_dir)
	
	# Index references for alignment
	print("Generating indexed DBs for alignment")
	index_ref(ref_dir, args.k_size, args.align_tool)
	# align and process
	print("Performing alignment")
	warn = align_baits(fasta_file, ref_dir, cov_dir, args.align_tool, args.input_db_file, threads, args.skip_plots)
	
	print("\nStarting statistical calculations")
	with open(os.path.join(cov_dir, "merged_coverage.txt"), 'w+') as out_merge:
		out_merge.write("Accession\tGenome Completeness\tAverage Depth of Coverage\tAverage Depth of Coverage Across Non-Zero\tBait Origin\tHit Origin(s)\n")
		for file in glob.glob(os.path.join(cov_dir, "*.txt")):
			if "merged_coverage.txt" == os.path.basename(file): continue
			for line in open(file):
				if line.startswith("Accession"): continue
				out_merge.write(line + "\n")
				if args.output_pre is False: os.remove(file)
	
	# ngsCAT coverage
	if args.run_ngscat is True:
		print("Running ngsCAT to calculate coverage statistics")
		ngscat(ref_dir, cov_dir, args.coveragethresholds, threads)
	else:
		print("Skipping ngsCAT analysis of coverage statistics.")
		
	if args.output_pre is False:
		clean(fasta_name, ref_dir)
	else:
		print("Intermediate files/directories will be kept:\n")
		print("Keeping directory %s" %(fasta_name+"_coverageref"))
		print("\nRemoving temporary files")
		for file in glob.glob(os.path.join(ref_dir, '*.bam')):
			filename = str(os.path.basename(file))
			if os.path.exists(os.path.join('/tmp', filename)): 
				os.remove(os.path.join('/tmp', filename))
	
	return warn


if __name__ == '__main__':
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	try:
		directory = script_path.rsplit("/", 1)[0]
	except IndexError:
		directory = os.getcwd() # script in current working directory
	warn = run()
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	#print("\nRemoved " + str(num_baits) + " bait(s).\n")
	print("\nStart: " + start)
	print("Finish: " + finish)
	if warn: print("WARNING: Some coverage information flagged missing! Check input databse file and re-run, if necessary!")
	exit()
