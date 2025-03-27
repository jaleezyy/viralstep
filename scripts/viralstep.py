#!/bin/env python

### MAIN BASE Viral Syndromic Target Enrichment Pipeline (ViralSTEP)
### classification - envoke classification branch to generate database file
###### generatedb (from FASTA file to database)
###### classify (using existing database file, generate database file for a given FASTA)
### vhost - envoke VHost-Classifier for given FASTA - filter for mammalian host
### baits - envoke bait design branch (flag for human whole genome capture)
###### all (from start to finish) - run default settings only
###### tiling ( + GC content)
###### tm
###### syndromic
###### crosshyb
###### redundacy
### stats - envoke statistics module for given FASTA (baits)
### evolution - envoke envolutionary bait design branch (in development)


import os, sys, shutil
import argparse
import subprocess
import tarfile, urllib.request

try:
	import scripts.classification as classify
	import scripts.vhost_filter as vhost
	import scripts.bait_coverage as coverage
	import scripts.syndromic_filter as syndromic
	import scripts.cross_hybridization as crosshyb
	import scripts.redundancy as redundant
	import scripts.human_only as human
	import scripts.random_syndromic_pull as subset
except ModuleNotFoundError: # if running viralstep.py directly
	import classification as classify
	import vhost_filter as vhost
	import bait_coverage as coverage
	import syndromic_filter as syndromic
	import cross_hybridization as crosshyb
	import redundancy as redundant
	import human_only as human
	import random_syndromic_pull as subset



class MainBase:
	def __init__(self):
		self.script_path = script_path = os.path.join(os.path.abspath(sys.argv[0]).rsplit("/",1)[0]) # script directory
		if not self.script_path.endswith("scripts"):
			if os.path.exists(os.path.join(self.script_path, "scripts")):
				self.script_path = os.path.join(self.script_path, "scripts") # executing from root of ViralSTEP
			else:
				exit("Cannot find ViralSTEP scripts directory!")
		self.HELP='''%(prog)s <command> [<args>]

Accepted commands (add -h/--help after command for more details):
	---------------------------------------------------------------------------------------
	Download Software, Packages, Public Databases
	---------------------------------------------------------------------------------------
	download         Download all required software, lanuage interpreters, and packages (i.e., R, Python, Ruby, Perl) 
	downloadvhost    Download VHost-Classifier and update corresponding Virus-Host DB file 
	downloadbaits    Download required software and packages for bait tiling and filtering (i.e., BaitsTools)
	downloadblastdb  Download (or update existing) BLASTn database (only nt)
	downloadstats    Download required software and packages for summary statistics and plots
	
	---------------------------------------------------------------------------------------
	Classifcation
	---------------------------------------------------------------------------------------
	updatetax        Update required taxonomy files for Taxonkit (will install to $HOME/.taxonkit)
	generatedb       From a given FASTA file, generate classification database file (.csv)
	classify         Using an existing database file (.csv), create a separate database file (.csv) for a given FASTA file

	---------------------------------------------------------------------------------------	
	VHost-Classifier (https://github.com/Kzra/VHost-Classifier; PMID: 30824917)
	---------------------------------------------------------------------------------------
	vhost_update     Update Virus-Host DB file (will default to VHost-Classifier location)
	vhost_filter     Invoke VHost-Classifier for a given FASTA file (filters for Mammalian)

	---------------------------------------------------------------------------------------	
	Bait Design
	---------------------------------------------------------------------------------------
	baits            Run default bait design protocol for a given FASTA file (SKIPS syndromicbaits and humanwgc)
	tilebaits        Run BaitsTools (tilebaits) for given input FASTA, bait length, offset, GC Content range
	tmbaits	         Run melting temperature filter for given FASTA and temperature range (inclusive)
	syndromicbaits   Run syndromic filter for given FASTA, classification database file, and desired syndrome
	crosshybaits     Run cross-hybridization filter for viral bait set
	redundbaits      Run filter to remove highly similar (sequence space) baits
	humanwgc         Run filter to remove non-human virus targetting baits (uses Virus-Host DB)
	normalize        Generate a randomized subset bait set with a maximum number of sequences per virus family (EXPERIMENTAL)
	
	---------------------------------------------------------------------------------------	
	Statistics
	---------------------------------------------------------------------------------------
	stats            Run all modules for bait statistics given a bait set FASTA file using default settings
	proportion       Breakdown bait set by virus family and Balitmore classifers (i.e., +ssRNA, dsDNA, etc.)
	physprop         Calculate GC Content and melting temperature for a given FASTA file with baits
	coverage         Calculate depth and breadth of coverage for a given bait set FASTA (optionally provided reference)
	plots            Generate summary virus family coverage plots
	
	---------------------------------------------------------------------------------------	
	General
	---------------------------------------------------------------------------------------
	help             Display this help page
	validate         Ensure all required packages and software are installed (will run through all functions with test data)
	
	'''
		parser = argparse.ArgumentParser(prog="viralstep", description='ViralSTEP: Viral Syndromic Target Enrichment Pipeline', usage=self.HELP)
		parser.add_argument('command', choices=['download', 'downloadvhost', 'downloadbaits', 'downloadblastdb', 'downloadstats',
												'updatetax', 'generatedb', 'classify', 
												'vhost_update', 'vhost_filter', 
												'baits', 'tilebaits', 'tmbaits', 'syndromicbaits', 'crosshybaits', 'redundbaits', 'humanwgc', 'normalize', 
												'stats', 'proportion', 'physprop', 'coverage', 'plots', 'help'],
												help='Subcommand to run')

		args=parser.parse_args(sys.argv[1:2])
		# if not hasattr(self, args.command):
			# exit("Error: Unrecognized command: {}".format(args.command))
		getattr(self, args.command)()

	def help(self):
		exit(self.HELP.replace("%s" %("%(prog)s"), "usage: viralstep"))
	
	def downloadstats(self):
		"""
		Download R packages mainly
		"""
	
	def downloadblastdb(self):
		"""
		Run update_blastdb.py. Will default --no-update on cross-hybridization filter.
		
		overwrite = recreate BLASTDB folder
		update = if found, check if any changes in database and update
		skip = if found, do not make any changes
		"""
		db_dir = os.path.join(self.script_path, 'blastdb')
		if os.path.isdir(db_dir):
			while True:
				answer = input("BLASTDB found at %s. Choose update action {'overwrite', 'update', 'skip'}: " %(db_dir))
				if answer.lower() in ('overwrite', 'update', 'skip'):
					action = answer.lower()
					break
				else:
					print("Invalid response!")
					continue
		if action == "overwrite":
			crosshyb.blastdb(True, False)
		if action == "update":
			crosshyb.blastdb(False, False)
		if action == "skip":
			crosshyb.blastdb(False, True)
		
	
	def downloadbaits(self):
		"""
		Pull all necessary git repos, packages
		"""
		pass

	def vhost_update(self):
		"""
		Consider importing function from human diagnostic filter
		"""
		vhost_dir = os.path.join(self.script_path, "VHost-Classifier")
		if not os.path.exists(vhost_dir):
			vhost_dir = os.getcwd()
		while True:
			answer = input("Select a directory to download database file to (leave blank for default: %s) :" %(vhost_dir))
			if os.path.isdir(answer): 
				dest = answer
				break
			else:
				print("Blank or invalid entry! Defaulting to VHost-Classifier location")
				dest = vhost_dir
				break
		human.download_db(dest)
			
	def downloadvhost(self):
		"""
		Pull from repo
		"""
		pass
	
	def download(self):
		"""
		Invoke individual download functions above
		"""
		pass
	
	def updatetax(self):
		"""
		Consider adding a means to update Taxonomy for KronaTools
		"""
		home_dir = os.getenv("HOME") # only works for UNIX-based devices
		#home_dir = os.getcwd()
		taxonkit_db = os.path.join(home_dir, '.taxonkit')
		req_files = ['names.dmp', 'nodes.dmp', 'delnodes.dmp', 'merged.dmp']
		
		if not os.path.exists(taxonkit_db):
			print("Creating required Taxonkit directory")
			os.mkdir(taxonkit_db)
		
		print("Downloading NCBI taxdump.tar.gz")
		filename, headers = urllib.request.urlretrieve("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", os.path.join(taxonkit_db, 'taxdump.tar.gz'))
		tax_tar = tarfile.open(os.path.join(taxonkit_db, 'taxdump.tar.gz'))
		print("Extracting required files")
		for file in req_files:
			tax_tar.extract(file, taxonkit_db)
		tax_tar.close()
		
		print("Cleaning up and setting environment")
		os.environ["TAXONKIT_DB"] = str(taxonkit_db)
		os.remove(os.path.join(taxonkit_db, 'taxdump.tar.gz'))
		exit("Required taxonomy database for Taxonkit updated!")
	
	def generatedb(self): 
		parser = self.generatedb_args()
		args = parser.parse_args(sys.argv[2:])
		self.generatedb_run(args)
	
	def generatedb_args(self):
		parser = argparse.ArgumentParser(prog="viralstep generatedb", description="Classification branch that take in a reference FASTA and generates a CSV database file to use for bait design.")
		parser.add_argument('-i', '--fasta', dest='input_fasta',required=True, help='Input FASTA file')
		parser.add_argument('-o', '--output', dest='output_dir', default='db_intermediate', help='Output directory for intermediate files', required=False)
		return parser

	def generatedb_run(self, args):
		script = os.path.join(self.script_path, 'generate_db.sh')
		if os.path.exists(script):
			subprocess.run(['bash', script, '-i', args.input_fasta, '-o', args.output_dir])
		else:
			print("Missing script 'generate_db.sh'")

	def classify(self):
		parser = self.classify_args()
		args = parser.parse_args(sys.argv[2:])
		self.classify_run(args)
	
	def classify_args(self):
		parser = classify.create_parser()
		return parser
		
	def classify_run(self, args):
		script = classify.__file__
		params = ['python3', script, '-i', args.input_fasta_file, \
									'-l', args.input_family_file, \
									'-o', args.output_file, \
									'-t', str(args.threads)]
		if args.receive_taxid is True: params.append('--taxid')
		if args.output_log is True: params.append('--log')
		subprocess.run(params)
	
	def vhost_filter(self):
		"""
		Currently divided into vhost_filter (main script) and vhost_parse (just to analysee VHost-Classifier data.
		TODO: Merge the functions into one file rather than have separate scripts
				- should just be able to take the function and re-work variables
		"""
		parser = self.vhost_filter_args()
		args = parser.parse_args(sys.argv[2:])
		self.vhost_filter_run(args)
		
	def vhost_filter_args(self):
		parser = vhost.create_parser()
		return parser
		
	def vhost_filter_run(self, args):
		script = vhost.__file__
		params = ['python3', script, '-i', args.input_fasta_file,\
									'-o', args.output_fasta_file,\
									'-c', args.input_db_csv
									]
		subprocess.run(params)
	
	def normalize(self):
		"""
		Experimental script creates randomized subset of a given bait set
		"""
		parser = self.normalize_args()
		args = parser.parse_args(sys.argv[2:])
		self.normalize_run(args)
	
	def normalize_args(self):
		parser = subset.create_parser()
		return parser
	
	def normalize_run(self, args):
		script = subset.__file__
		
		params = ['python3', script, '-i', args.input_fasta_file,\
									'-db', args.db_file,\
									'-o', args.output_fasta_file,\
									'-s', args.pull_length
									]
		if args.subset_db is True: params.append('--generatedb')
		subprocess.run(params)
	
	def humanwgc(self):
		"""
		Check if VirusHostDB file exists in VHost-Classifier by default, otherwise run vhost_update
		(Not to be included in 'baits' function below)
		"""
		parser = self.humanwgc_args()
		args = parser.parse_args(sys.argv[2:])
		self.humanwgc_run(args)
		
	def humanwgc_args(self):
		parser = human.create_parser()
		return parser
		
	def humanwgc_run(self, args):
		script = human.__file__
		
		### If --database is None: check VHost-Classifier location, otherwise specify where or run vhost_update
		if args.virus2hostdb is None:
			vhost_db = os.path.join(self.script_path, "VHost-Classifier", "virushostdb.tsv")
		else:
			vhost_db = args.virus2hostdb
		if not os.path.exists(vhost_db): exit("Cannot find Virus-Host DB file. Either specify file location or run 'viralstep vhost_update'")
		
		params = ['python3', script, '-i', args.input_fasta_file,\
									'-db', vhost_db,\
									'-l', args.taxid_list,\
									'-o', args.output_fasta_file
									]
		if args.output_pre is True: params.append('--keep-intermediate')
		subprocess.run(params)
	
	def redundbaits(self):
		"""
		binb4greedy AKA redundancy filter
		"""
		parser = self.redundbaits_args()
		args = parser.parse_args(sys.argv[2:])
		self.redundbaits_run(args)
		
	def redundbaits_args(self):
		parser = redundant.create_parser()
		return parser
	
	def redundbaits_run(self, args):
		script = redundant.__file__
		params = ['python3', script, '-i', args.input_fasta_file,\
									'-o', args.output_fasta_file,\
									'-t', str(args.num_threads)
									]
		if args.align_length is not None: params.append('-l %s' %(args.align_length))
		if args.align_sim is not None: params.append('-s %s' %(args.align_sim))
		if args.sequence_ident is not None: params.append('-m %s' %(args.sequence_ident))
		if args.output_pre is True: params.append("--keep-intermediate")
		subprocess.run(params)
	
	def crosshybaits(self):
		"""
		If BLASTDB found, use it. Otherwise run downloadblastdb
		"""
		parser = self.crosshybaits_args()
		args = parser.parse_args(sys.argv[2:])
		self.crosshybaits_run(args)
		
	def crosshybaits_args(self):
		parser = crosshyb.create_parser()
		return parser
		
	def crosshybaits_run(self, args):
		script = crosshyb.__file__
		params = ['python3', script, '-i', args.input_fasta_file,\
									'-o', args.output_fasta_file,\
									'-t', str(args.num_threads),\
									'-k', args.keep_term
									]
		if args.input_blast_file is not None: params.append('-b %s' %(args.input_blast_file))
		if args.validate_baits is True: params.append('--validate')
		if args.only_validate is True: params.append('--validate-only')
		if args.output_pre is True: params.append('--keep-intermediate')
		if args.strict_filter is True: params.append('--strict')
		subprocess.run(params)
	
	def syndromicbaits(self):
		"""
		Do not include in 'baits' function below. 
		TODO: Add ability to submit comma-separated list of virus families for filtration
		TODO: Consider re-working handling of NCBI-styled headers
		"""
		parser = self.syndromicbaits_args()
		args = parser.parse_args(sys.argv[2:])
		self.syndromicbaits_run(args)
	
	def syndromicbaits_args(self):
		parser = syndromic.create_parser()
		return parser
		
	def syndromicbaits_run(self, args):
		script = syndromic.__file__
		params = ['python3', script, '-i', args.input_fasta_file,\
									'-l', args.input_family_file,\
									'-o', args.output_fasta_file,\
									'-f', args.filter_name,\
									'-t', args.threads
									]
		if args.not_filtered is True: params.append('--rejected')
		if args.receive_taxid is True: params.append('--taxid')
		if args.receive_lineage is True: params.append('--lineage')
		if args.output_log is True: params.append('--log')
		subprocess.run(params)
	
	def tmbaits(self):
		"""
		Melting temperature filter runs here. Needs to connect with melt_parse_baits.sh.
		"""
		parser = self.tmbaits_args()
		args = parser.parse_args(sys.argv[2:])
		self.tmbaits_run(args)
	
	def tmbaits_args(self):
		parser = argparse.ArgumentParser(prog="viralstep tmbaits", description="Melting temperature filter for a given bait set.")
		parser.add_argument('-i', '--fasta', dest='input_fasta',required=True, help='Input FASTA file')
		parser.add_argument('-min', '--mintemp', dest='min_temp', default=65, type=float, help='Minimum melting temperature in degrees Celcius (inclusive) (Default = 65.0)')
		parser.add_argument('-max', '--maxtemp', dest='max_temp', default=95, type=float, help='Maximum melting temperature in degrees Celcius (inclusive) (Default = 95.0)')
		parser.add_argument('-tm', '--tmlist', dest='tm_list', default=None, help="Existing melting temperature list. Can be used to skip calculations and go straight to filter.")
		return parser
		
	def tmbaits_run(self, args): ### TODO: Adjust for proper import
		script = os.path.join(self.script_path, 'melting_temperature.py')
		### Check if required script is present
		script_melt_parser = os.path.join(self.script_path, 'melt_parse_baits.sh')
		if ((not os.path.exists(script_melt_parser)) and (args.tm_list is None)): 
			exit("Missing melt_parse_baits.sh to calculate melting temperatures. Please check scripts or supply melting temperature list --tmlist")
		
		### Create symlink if --tmlist provided
		clean = False
		if ((args.tm_list is not None) and (os.path.exists(args.tm_list))):
			print("Creating symlink to melting temperature list")
			dest_file = os.path.join(os.path.dirname(args.input_fasta), "Tm_{}.txt".format(os.path.basename(args.input_fasta)))
			os.symlink(args.tm_list, dest_file)
			clean = True
		
		params = ['python3', script, args.input_fasta, str(args.min_temp), str(args.max_temp)]
		subprocess.run(params)
		if clean == True: os.remove(dest_file)
	
	def tilebaits(self):
		"""
		Requires BaitsTools (v1.6.7+). GC Content filter here. 
		"""
		parser = self.tilebaits_args()
		args = parser.parse_args(sys.argv[2:])
		self.tilebaits_run(args)
	
	def tilebaits_args(self): ### Match with BaitsTools arguments 
		USAGE="""%(prog)s [<args>]

Default tiling (flags highlighted will be preset but others such as --threads and output flags may still be used):

	    --default				Parameters used: --length 80 --offset 20 --complete --noNs --mingc 25.0 --maxgc 55.0
	    --default-ncbi			Same as --default but with added -D (or --ncbi) for assumed NCBI-styled descriptions

tilebaits-specific options:

	-i, --input [FILE]			Input FASTA/FASTQ file (REQUIRED)
	-L, --length [VALUE]			Requested bait length (Default = 120)
	-O, --offset [VALUE]			Base pair offset between tiled baits (Default = 60)

Bait filtration options:

	-c, --complete				Require baits be full length
	-N, --noNs				Exclude bait sequences with Ns
	-n, --mingc [VALUE]			Minimum GC content in %% (Default = 30.0)
	-x, --maxgc [VALUE]			Maximum GC content in %% (Default = 50.0)

Common options:	

	-o, --outprefix [VALUE]			Output file prefix (Default = out)
	-Z, --outdir [VALUE]			Output directory path (Default is ./)
	-l, --log				Output detailed log
	-B, --bed				Output BED file for the baits
	-E, --rbed				Output BED file for the baits relative to extracted sequences
	    --shuffle				Shuffle baits to compensate for extending beyond contig ends
	-D, --ncbi				FASTA/FASTQ file headers have NCBI-style descriptions
	-C, --collapse				Collapse ambiguities to a single nucleotide
	-Y, --rna				Output baits as RNA rather than DNA
	-R, --rc				Output reverse-complemented baits
	-G, --gaps [VALUE]			Strategy to handle sequence gaps (-) (include, exclude, or extend) (Default = include)
	-5, --5prime [VALUE]			Sequence to addend to 5' end of baits
	-3, --3prime [VALUE]			Sequence to addend to 3' end of baits
	-X, --threads [VALUE]			Number of threads (Default = 1)
	    --gzip				Gzip output files
"""
		parser = argparse.ArgumentParser(prog="viralstep tilebaits", description="Envoke BaitsTools (v1.6.7+) to tile baits from a given FASTA/FASTQ.", usage=USAGE)
		### default settings
		parser.add_argument('--default', dest='default', action='store_true', help="Run default parameters (see usage for details)")
		parser.add_argument('--default-ncbi', dest='defaultncbi', action='store_true', help="Run default parameters assuming NCBI-style descriptions (see usage for details)")
		### Input
		parser.add_argument('-i', '--input', dest='input', required=True, help='Input FASTA/FASTQ file')
		parser.add_argument('-L', '--length', dest='length', default=120, type=int, help='Requested bait length (Default = 120)')
		parser.add_argument('-O', '--offset', dest='offset', default=60, type=int, help='Base pair offset between tiled baits (Default = 60)')
		### Filters
		parser.add_argument('-c', '--complete', dest='complete', action='store_true', help='Require baits be full length')
		parser.add_argument('-N', '--noNs', dest='noNs', action='store_true', help='Exclude bait sequences with Ns')
		parser.add_argument('-n', '--mingc', dest='minGC', default=30.0, type=float, help='Minimum GC content in %% (Default = 30.0)')
		parser.add_argument('-x', '--maxgc', dest='maxGC', default=50.0, type=float, help='Maximum GC content in %% (Default = 50.0)')
		### Output
		parser.add_argument('-o', '--outprefix', dest='outprefix', default='out', help='Output file prefix (Default = out)')
		parser.add_argument('-Z', '--outdir', dest='outdir', default=os.getcwd(), help='Output directory path (Default is %s)' %(os.getcwd()))
		### Settings
		parser.add_argument('-l', '--log', dest='outlog', action='store_true', help="Output detailed BaitsTools log")
		parser.add_argument('-B', '--bed', dest='outbed', action='store_true', help="Output BED file for the baits")
		parser.add_argument('-E', '--rbed', dest='outrbed', action='store_true', help="Output BED file for the baits relative to extracted sequences")
		parser.add_argument('-D', '--ncbi', dest='ncbi_headers', action='store_true', help="FASTA/FASTQ file headers have NCBI-style descriptions")
		parser.add_argument('-C', '--collapse', dest='collapse', action='store_true', help="Collapse ambiguities to a single nucleotide")
		parser.add_argument('-Y', '--rna', dest='outrna', action='store_true', help="Output baits as RNA rather than DNA")
		parser.add_argument('-R', '--rc', dest='rev_complement', action='store_true', help="Output reverse-complemented baits")
		parser.add_argument('-G', '--gaps', dest='handle_gaps', default='include', choices=['include', 'exclude', 'extend'], help = "Strategy to handle sequence gaps (-) (include, exclude, or extend) (Default = include)")
		parser.add_argument('-5', '--5prime', dest='fiveprime_seq', help="Sequence to addend to 5' end of baits")
		parser.add_argument('-3', '--3prime', dest='threeprime_seq', help="Sequence to addend to 3' end of baits")
		parser.add_argument('-X', '--threads', dest='threads', default=1, type=int, help="Number of threads (Default = 1)")
		parser.add_argument('--gzip', dest='outgzip', action='store_true', help="Gzip output files")
		return parser
		
	def tilebaits_run(self, args):
		### Test BaitsTools prior to running
		baitstools = os.path.join(self.script_path, 'BaitsTools', 'baitstools.rb')
		try: #subprocess.call can be used instead, remove check=True
			subprocess.run(['ruby', baitstools, '-v'], stdout=subprocess.DEVNULL, check=True)
			subprocess.run(['ruby', baitstools, 'tilebaits', '-v'], stdout=subprocess.DEVNULL, check=True)
		except CalledProcessError:
			exit("Something went wrong when executing BaitsTools. Ensure the software is located %s" %(os.path.dirname(baitstools)))
		
		params = ['ruby', baitstools, 'tilebaits', '-i', args.input,\
													'-L', str(args.length),\
													'-O', str(args.offset),\
													'--mingc', str(args.minGC),\
													'--maxgc', str(args.maxGC),\
													'--outprefix', args.outprefix,\
													'--outdir', args.outdir,\
													'--gaps', args.handle_gaps,\
													'--threads', str(args.threads)
													]
		### If default, override args
		if ((args.default is True) and (args.defaultncbi is False)):
			params = ['ruby', baitstools, 'tilebaits', '-i', args.input,\
													'-L', '80',\
													'-O', '20',\
													'--mingc', '25.0',\
													'--maxgc', '55.0',\
													'--outprefix', args.outprefix,\
													'--outdir', args.outdir,\
													'--gaps', args.handle_gaps,\
													'--threads', str(args.threads),\
													'--complete', '--noNs'
													]
		elif ((args.default is False) and (args.defaultncbi is True)):
			params = ['ruby', baitstools, 'tilebaits', '-i', args.input,\
													'-L', '80',\
													'-O', '20',\
													'--mingc', '25.0',\
													'--maxgc', '55.0',\
													'--outprefix', args.outprefix,\
													'--outdir', args.outdir,\
													'--gaps', args.handle_gaps,\
													'--threads', str(args.threads),\
													'--ncbi', '--complete', '--noNs'
													]
		elif ((args.default is True) and (args.defaultncbi is True)):
			print("Both default options present! Running --default")
			params = ['ruby', baitstools, 'tilebaits', '-i', args.input,\
													'-L', '80',\
													'-O', '20',\
													'--mingc', '25.0',\
													'--maxgc', '55.0',\
													'--outprefix', args.outprefix,\
													'--outdir', args.outdir,\
													'--gaps', args.handle_gaps,\
													'--threads', str(args.threads),\
													'--complete', '--noNs'
													]
		else: # both --default and --defaultncbi are False
			pass
		if ((args.complete is True) and ('--complete' not in params)): params.append('--complete')
		if ((args.noNs is True) and ('--noNs' not in params)): params.append('--noNs')
		if args.outlog is True: params.append('--log')
		if args.outbed is True: params.append('--bed')
		if args.outrbed is True: params.append('--rbed')
		if ((args.ncbi_headers is True) and ('--ncbi' not in params)): params.append('--ncbi')
		if args.collapse is True: params.append('--collapse')
		if args.outrna is True: params.append('--rna')
		if args.rev_complement is True: params.append('--rc')
		if args.fiveprime_seq is not None: params.extend(['--5prime', args.fiveprime_seq])
		if args.threeprime_seq is not None: params.extend(['--3prime', args.threeprime_seq])
		if args.outgzip is True: params.append('--gzip')

		subprocess.run(params)
	
	def baits(self): ### May remove
		"""
		Run all bait steps (above) with default settings
		"""
		pass
	
	def plots(self):
		"""
		Will keep zeroes by default"
		"""
		parser = self.plots_args()
		args = parser.parse_args(sys.argv[2:])
		self.plots_run(args)
		
	def plots_args(self):
		USAGE="""%(prog)s [<args>]

Syndromes accepted include: 
	all (default)
	respir (or respiratory)
	blood
	gastro (or gastrointestinal)
	human (for all human-relevant virus families)
	csf (meaning cerebral spine fluid)
	urine
	sjf (meaning skin or joint fluid)
	special (referring to special pathogens, zoonotics, vectorborne)

Use -h/--help for more details.
	"""
		parser = argparse.ArgumentParser(prog="viralstep plots", description="Generate summary plots.", usage=USAGE)
		parser.add_argument('-o', '--output', dest='output', required=True, help='Output directory name (REQUIRED)')
		parser.add_argument('-p', '--phys-prop', dest='phys_prop', default="NULL", help='File containing GC content and melting temperature (in that order) for a bait set in TSV format')
		parser.add_argument('-f', '--proportion', dest='proportion', default="NULL", help='File containing the number (and percentage) of baits per virus family in TSV format (family, number of baits, percentage of baits)')
		parser.add_argument('-b', '--baltimore', dest='baltimore', default="NULL", help='File containing the number (and percentage) of baits per balitmore classifier (ex. +ssRNA, dsDNA, etc.) in TSV format (classifier, number of baits, percentage of baits)')
		parser.add_argument('-c', '--coverage', dest='coverage', default="NULL", help='File containing coverage information in TSV format (accession, genome completeness, average depth of coverage)')
		parser.add_argument('-s', '--syndrome', dest='syndrome', default='all', help='Examine coverage for specific groups of virus families (or syndromes). Input can be a comma-separated list of virus families or one of the accepted values (see usage above).')
		parser.add_argument('--remove-zero', dest='remove_zero', action='store_true', help="Remove 0 coverage information from plots and tables (only applies if '--coverage' is provided)")
		parser.add_argument("--non-zero", dest='non_zero', action="store_true", help="Assess (depth of) coverage across non-zero regions only. Useful for gapped coverages where breadth of coverage (completeness) is expected to be low")
		return parser
		
	def plots_run(self, args):
		script = os.path.join(self.script_path, 'plots.R')
		params = ['Rscript', script, '--output', args.output, \
									'--phys-prop', args.phys_prop, \
									'--proportion', args.proportion, \
									'--baltimore', args.baltimore, \
									'--coverage', args.coverage, \
									'--syndrome', args.syndrome]
		if args.remove_zero: params.append('--remove-zero')
		if args.non_zero: params.append('--non-zero')
		if os.path.exists(script):
			subprocess.run(params)
		else:
			print("Missing script 'plots.R'")
	
	def coverage(self):
		"""
		Will use Bowtie2 alignment by default"
		"""
		parser = self.coverage_args()
		args = parser.parse_args(sys.argv[2:])
		self.coverage_run(args)
		
	def coverage_args(self):
		parser = coverage.create_parser()
		return parser
		
	def coverage_run(self, args):
		script = coverage.__file__
		params = ['python3', script, '-i', args.input_fasta_file, \
									'-d', args.input_db_file, \
									'-r', args.manual_ref, \
									'-a', args.align_tool, \
									'-t', str(args.num_threads), \
									'--coveragethrs', args.coveragethresholds, \
									'-k', args.k_size]
		if args.run_ngscat is True: params.append('--ngscat')
		if args.output_pre is True: params.append('--keep-intermediate')
		if args.skip_plots is True: params.append('--skip-plots')
		subprocess.run(params)
	
	def physprop(self):
		"""
		Calculate GC Content and melting temperature (can take up excess space due to melt.pl)
		"""
		parser = self.physprop_args()
		args = parser.parse_args(sys.argv[2:])
		self.physprop_run(args)
	
	def physprop_args(self):
		parser = argparse.ArgumentParser(prog="viralstep physprop", description="Calculate GC content and melting temperatures for given bait set in FASTA format.")
		parser.add_argument('-i', '--fasta', dest='input_fasta', required=True, help='Input FASTA file')
		parser.add_argument('-o', '--output', dest='output_file', default='phys_prop.txt', help='Output file. Default: phys_prop.txt')
		parser.add_argument('--keep', dest='keep_intermediate', action='store_true', help='Keep intermediate GC content and melting temperature data')
		return parser

	def physprop_run(self, args): ### TODO: Adjust for proper import
		script = os.path.join(self.script_path, 'phys_prop.py')
		
		### Check if required scripts (GCcontent.py and melt_parse_baits.sh are present)
		
		
		params = ['python3', script, args.input_fasta, args.output_file]
		params.append(str(args.keep_intermediate)) # can be True or False
		if os.path.exists(script):
			subprocess.run(params)
		else:
			print("Missing script 'phys_prop.py'")

	def proportion(self):
		"""
		Output files by default
		"""
		parser = self.proportion_args()
		args = parser.parse_args(sys.argv[2:])
		self.proportion_run(args)
	
	def proportion_args(self):
		parser = argparse.ArgumentParser(prog="viralstep proportion", description="Calculate proportions of virus families and baltimore classifiers from provided CSV database file.")
		parser.add_argument('-i', '--inputdb', dest='input_db', required=True, help='Input classification CSV file')
		parser.add_argument('-o', '--output', dest='output_prefix', default=None, help="Output file prefix. If not provided, default will output to stdout")
		return parser
		
	def proportion_run(self, args): ### TODO: Adjust for proper import
		script = os.path.join(self.script_path, 'proportion.py')
		if args.output_prefix is not None: 
			params = ['python3', script, args.input_db, args.output_prefix]
		else:
			params = ['python3', script, args.input_db]
		if os.path.exists(script):
			subprocess.run(params)
		else:
			print("Missing script 'proportion.py'")
	
	def stats(self): ### May remove
		"""
		Run all statistical modules (above) with default settings
		"""
		exit("Run individually for now...")
	

if __name__ == '__main__':
	MainBase()
