#!/usr/bin/env python

### Goal is to provide FASTA input from clustering with an incomplete header and update said header with taxonomy information
### Output header would be <accessionID>_<virus_Family>
### Uses ViralSTEP DB to convert accID --> virus family
### Top hit of alignment via KMA >> Bowtie2 used 

from dev_scripts import alignment
from Bio import SeqIO
import argparse
import os, sys
from shutil import rmtree

	
def family_hash(db):
	accID2family = {}
	with open(db) as f:
		for line in f:
			col = [item.strip() for item in line.split(",")]
			accID2family[col[0]] = col[-1] # first item contains accID, last contains virus family
		
	return accID2family
	
def cleanup(tmp):
	tmp_no_dup = set(tmp)
	for file in tmp_no_dup:
		assert os.path.exists(file)
		if os.path.isdir(file):
			rmtree(file)
		elif os.path.isfile(file):
			os.remove(file)
		else:
			print(f"Could not remove {file}")
	
def run(args):
	### Parameters
	fasta_file = os.path.abspath(args.input_fasta_file)
	fasta_name = os.path.basename(fasta_file)
	ref_file = os.path.abspath(args.manual_ref)
	ref_name = os.path.basename(ref_file)
	threads = int(args.threads)
	kmer = int(args.kmer)
	
	
	tmp_files = []
	
	### determine output
	if args.output_fasta_file is None:
		final_output = f"corrected_{fasta_name}"
	else:
		final_output = args.output_fasta_file
		
	if os.path.exists(final_output) and not args.overwrite:
		print(f"{final_output} cannot already exist! Use --overwrite or rename output!")
		sys.exit(1)
	elif os.path.exists(final_output) and args.overwrite:
		os.remove(final_output)
	else:
		pass # no action needed
		
	### Generate conversion hash for accID --> virus family
	if args.input_db_file is not None and os.path.exists(args.input_db_file):
		print("ViralSTEP DB found!")
		db_file = os.path.abspath(args.input_db_file)
		db_name = os.path.basename(db_file)
		accID2family = family_hash(db_file)
	else:
		accID2family = None
		# will use query and ref names in alignment if not provided
	
	### Create indexes for reference sequence(s)
	if not os.path.exists(ref_file):
		print(f"MISSING: {ref_name} cannot be found!")
		sys.exit(1)
	else:
		print("Input reference FASTA found!")
		ref_index, ref_index_out = alignment.index_ref(ref_file, os.getcwd(), kmer)
		tmp_files.extend(ref_index_out)
		
	### Run KMA to produce alignment
	if not os.path.exists(fasta_file):
		print(f"MISSING: {fasta_name} cannot be found!")
	else:
		align_intermediates = {} # stored as bam_file: consensus_file
		align_seq = {} # sotred as bam_file: seq
		total_seq = alignment.count_seq(fasta_file)
		print(f"{total_seq} sequences found!")
		
		# split FASTA --> alignment is one query --> many references
		for seq_record, count in zip(SeqIO.parse(fasta_file, 'fasta'), range(1, total_seq+1)):
			id = str(seq_record.id) # shortened description
			header = seq_record.description
			seq = str(seq_record.seq).upper()
			
			tmp_fasta = os.path.join(os.getcwd(), f"{id}_{count}.fasta.tmp")
			
			# output tmp FASTA
			with open(tmp_fasta, 'w+') as outtmp:
				outtmp.write(f">{header}\n{seq}")
			
			tmp_files.append(tmp_fasta)
			# verify single sequence output
			assert alignment.count_seq(tmp_fasta) == 1
			
			# run KMA
			bam_file, consensus_file = alignment.align_seqs(tmp_fasta, os.getcwd(), ref_index, threads)
			
			# add output to collective dict and note sequence per bam
			align_intermediates[bam_file] = consensus_file
			align_seq[bam_file] = (header, seq)
		
			# note tmp files (to be deleted)
			tmp_files.extend([bam for bam in align_intermediates])
			tmp_files.extend([align_intermediates[con] for con in align_intermediates])
		
		# generate collective output file per sequence after above loop
		with open(final_output, 'w+') as out_fasta:
			for bam in align_intermediates:
				# check consensus file (.sam.fsa), if 1 sequence found, assign
				if alignment.count_seq(align_intermediates[bam]) == 1:
					record = SeqIO.read(align_intermediates[bam], 'fasta')
					accession = str(record.id)
					
				else:
					# check alignment (.bam) for top/only hit?
					# presumably there will be only one hit (to the intended reference)
					# otherwise check all hits and figure out top hit
					# perhaps count hit per reference
					accession = alignment.process_align(bam)
				
				# determine virus family
				if accID2family is not None:
						try:
							family = accID2family[accession]
						except KeyError:
							family = "Unknown"
				else:
					family = "Unknown"
					
				# generate updated header
				update_header = f"{accession}_{family}"
				update_seq = align_seq[bam][1]
				
				# check if we're keeping old header, append
				# replace existing "_" with "-" for compatibility
				if args.keep_old_headers:
					legacy_header = str(align_seq[bam][0]).replace("_", "-")
					update_header = f"{update_header}_{legacy_header}"
				
				# final output for sequence
				out_fasta.write(f">{update_header}\n{update_seq}\n")
			
		
	if not args.keep_intermediates:
		cleanup(tmp_files)
	else:
		print("NOTE: All intermediate files have not been deleted!")

def create_parser():
	parser = argparse.ArgumentParser(description='Perform alignment per sequence in a given cluster to reference sequences and determine likely taxonomy. Headers will be updated through this script.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="Input FASTA file with sequences needing taxonomic identification")
	parser.add_argument('-o', '--output', dest='output_fasta_file', required=False, default=None, help='Output FASTA file. By default an output filename will be generated based on the input file')
	parser.add_argument('-d', '--database', dest='input_db_file', required=False, default=None, help="ViralSTEP classificaton database file.")
	parser.add_argument('-r', '--reference', dest="manual_ref", required=True, help="Single FASTA file containing all desired reference sequences.")
	#parser.add_argument('-a', '--aligner', dest="align_tool", required=False, default="kma", help="Specify aligner to use between KMA/Bowtie2. Default=KMA")
	parser.add_argument('-k', '--kmer', dest='kmer', type=int, default=16, help='Optional argument for KMA aligner. Will index and align accordingly. Default k-mer size of 16 will be used if not specified.')
	parser.add_argument('--overwrite', dest='overwrite', action='store_true', help="Overwrite any existing intermediate or temporary files/directories")
	parser.add_argument('--retain-headers', dest='keep_old_headers', action='store_true', help='If toggled, then old header names will be maintained in the final corrected output. Underscore characters (_) will be replaced with hyphens (-).')
	parser.add_argument('--keep-intermediates', action='store_true', help="For debugging purposes, keep all intermediate files.")
	parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads. Default is 1.")
	return parser
	
if __name__ == '__main__':
	script_path = sys.argv[0]
	start_dir = os.getcwd()
	tmp_files = []
	
	try:
		script_dir = script_path.rsplit("/", 1)[0]
	except IndexError:
		script_dir = os.getcwd() # script in current working directory
		
	parser = create_parser()
	args = parser.parse_args()
	
	global threads
	threads = args.threads
	global kmer_size
	kmer_size = args.kmer
	
	run(args)
	sys.exit(0)