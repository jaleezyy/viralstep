#!/bin/env python

### Goal is to run the necessary programs and scripts within ViralSTEP to generate a table of GC Content and Melting Temperature

### Will generate .tsv (.txt) file using GCcontent.py
### Will run melt_parse_baits.sh to generate .tsv (.txt) of melting temperature

### Final step is to merge column values between two tab-delimited files
### Merge will occur in the order of provided files (i.e. columns of <input_file_2> will be appended to columns AFTER <input_file_1>)

import sys, os
import subprocess

HELP = """
Arguments are as followed:
phys_prop.py <fasta_file> <output_file> <keep_intermediate [True/*False*]>
		
Note: keep_intermediate == False by default. If True, script will only keep GC Content and Melting Temperature script output.
"""

try:
	fasta_input = sys.argv[1] 
	output_file = str(sys.argv[2])
	keep = eval(sys.argv[3])
except (IndexError, NameError): # first two arguments required
	try:
		fasta_input = sys.argv[1] 
		output_file = str(sys.argv[2])
		keep = False # accounts for NameError
	except (IndexError, NameError):
		print(HELP)
		exit()

def run_gc(script, fasta, name):
	print("\nCalculating GC Content")
	with open("GC_%s.txt" %(name), "w+") as gc_out:
		subprocess.run([script, fasta], stdout=gc_out, universal_newlines = True)
	
	return "GC_%s.txt" %(name)
	
def run_tm(script, fasta, name):
	print("Running melt.pl to calculate melting temperatures")
	with open("Tm_%s.txt" %(name), "w+") as tm_out:
		subprocess.run([script, "-s", fasta, "-o", "Tm_%s.txt" %(name)], stdout=tm_out, universal_newlines = True)
	
	return "Tm_%s.txt" %(name)


def output(gc, tm, output_file):
	print("Merging physical properties data")
	with open(gc) as f, open(tm) as f2:
		assert(len(f.readlines()) == len(f2.readlines()))
	with open(gc) as f, open(tm) as f2, open(output_file, "w+") as out:
		for line1,line2 in zip(f,f2):
			out.write(line1.strip("\n") + "\t" + line2.strip("\n") + "\n")

def clean(tmp_files, fastafile, keep):
	print("Cleaning up")
	if keep: print("Keeping intermediate GC and melting temperature data files")
	to_keep = ["Tm_%s.txt" %(fastafile), "GC_%s.txt" %(fastafile)]
	for file in tmp_files:
		if (keep is True) and (file in to_keep):
			continue
		else:
			os.remove(file)

### Startup parameters
tmp_files = []
gc_output = None
tm_output = None

try:
	fastafile = fasta_input.rsplit("/",1)[1]
	fastapath = fasta_input.rsplit("/",1)[0]
except IndexError: # current working directory
	fastafile = fasta_input
	fastapath = os.getcwd()

### Verify scripts present 
### TODO: Add download for equivalent script if script not found, add to tmp_files
gc_path = os.path.join(os.path.abspath(sys.argv[0]).rsplit("/",1)[0], "GCcontent.py")
tm_path = os.path.join(os.path.abspath(sys.argv[0]).rsplit("/",1)[0], "melt_parse_baits.sh")

### GC Content
if os.path.exists(gc_path):
	if os.path.exists("%s/GC_%s.txt" %(fastapath, fastafile)) :
		print("GC Content data exists for original FASTA file, skipping GCcontent.py" + "\n")
		gc_output = os.path.join(fastapath,"GC_%s.txt" %(fastafile))
	elif os.path.exists("GC_%s.txt" %(fastafile)):
		print("GC Content data exists for original FASTA file, skipping GCcontent.py" + "\n")
		gc_output = "GC_%s.txt" %(fastafile)
	else:
		gc_output = run_gc(gc_path, fasta_input, fastafile)
		tmp_files.append("GC_%s.txt" %(fastafile))
else:
	print("Missing GCcontent.py")
	### TODO: Add download link

### Melting Temperature
if os.path.exists(tm_path):
	if os.path.exists("%s/Tm_%s.txt" %(fastapath, fastafile)):
		print("Melting temperatre data exists for original FASTA file, skipping melt_parsh.sh" + "\n")
		tm_output = os.path.join(fastapath,"Tm_%s.txt" %(fastafile))
	elif os.path.exists("Tm_%s.txt" %(fastafile)):
		print("Melting temperatre data exists for original FASTA file, skipping melt_parsh.sh" + "\n")
		tm_output = "Tm_%s.txt" %(fastafile)
	else:
		tm_output = run_tm(tm_path, fasta_input, fastafile) # add files to tmp_files
		tmp_files.extend(["Tm_%s.txt" %(fastafile), "%s.symlink.65.ext" %(fastafile), "%s.symlink.65.plot" %(fastafile), \
							"%s.symlink.ct" %(fastafile), "%s.symlink.dG" %(fastafile), "%s.symlink.run" %(fastafile)])
		
else:
	print("Missing melt_parse_baits.sh")
	### TODO: Add download link

### Merge into output
if (gc_output is None) and (tm_output is None):
	print("Something went wrong! No existing GC or Melting Temperature data or scripts to calculate physical properties found!")
	exit()
elif (gc_output is not None) and (tm_output is None):
	print("Only GC Content data found or calculated: %s" %(gc_output))
elif (gc_output is None) and (tm_output is not None):
	print("Only Melting Temperature data found or calculated: %s" %(tm_output))
else:
	try:
		output(gc_output, tm_output, output_file)
	except AssertionError: # uneven files, keep intermediates (override) and try again
		print("\nGC content and Melting Temperature found or calculated, but is uneven. No merged file output will be provided.")
		print("Check intermediate files and try again!\n")
		clean(tmp_files, fastafile, True)
		exit()

### Clean up
if len(tmp_files) > 0: clean(tmp_files, fastafile, keep)

print("Done\n") 
exit()
			
