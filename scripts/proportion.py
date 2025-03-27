#!/bin/env python

### Goal is to take classification files from syndromic filter (or standalone classification.py)
### Input file contains <taxID>,<virus family>
### And calculate proportion for each virus family
### Mitigates the need for immediate R analysis for plots
### Gives rough numbers

### TODO: Add matplotlib options for virus family and baltimore breakdown

import sys

HELP="""
Arguments are as followed:
	proportion.py <input_file> <output_prefix>
	
	Input file should be a CSV with the following columns: AccessionID, TaxID, Baltimore classification, Virus family.
	If no filename given, output will only go to stdout.
"""

try:
	if "help" in sys.argv: exit(HELP)
	input = sys.argv[1]
	prefix = sys.argv[2]
	write = True
	#index_number = sys.argv[2]
	#output = sys.argv[2]
except IndexError:
	try:
		if help in sys.argv: exit(HELP)
		input = sys.argv[1]
		prefix = None
		write = False
	except IndexError:
		exit(HELP)
	
def count_line(input):
	num_line = len([1 for line in open(input)])
	return num_line
	
num_line = count_line(input)
family = {}
baltimore = {}
broad_balt = {}

with open(input, 'r') as f:
	for i in f:
		v_family = str(i).split(",")[-1]
		balt = str(i).split(",")[-2]
		balt_broad = balt.strip("+").strip("-")
		if v_family not in family:
			family[v_family] = 1
		else:
			family[v_family] = family[v_family] + 1
		if balt not in baltimore:
			baltimore[balt] = 1
		else:
			baltimore[balt] = baltimore[balt] + 1
		if balt_broad not in broad_balt:
			broad_balt[balt_broad] = 1
		else:
			broad_balt[balt_broad] = broad_balt[balt_broad] + 1

# virus family breakdown
if write:
	virus_family = open("%s_virus_family.tsv" %(prefix), 'w+')
	virus_baltimore = open("%s_baltimore.tsv" %(prefix), 'w+')
	virus_broad = open("%s_broad_baltimore.tsv" %(prefix), 'w+')

for f in family:
	print(str(f).split("\n")[0] + "\t" + \
	str(float(family[f]/num_line)*100).split("\n")[0] + "%" + "\t" + \
	str(int(family[f])) + " matches")
	if write:
		virus_family.write(str(f).split("\n")[0] + "\t" + \
		str(float(family[f]/num_line)*100).split("\n")[0] + "%" + "\t" + \
		str(int(family[f])) + " matches\n")
	
print("\n")

for b in baltimore:
	print(str(b).split("\n")[0] + "\t" + \
	str(float(baltimore[b]/num_line)*100).split("\n")[0] + "%" + "\t" + \
	str(int(baltimore[b])) + " matches")
	if write:
		virus_baltimore.write(str(b).split("\n")[0] + "\t" + \
		str(float(baltimore[b]/num_line)*100).split("\n")[0] + "%" + "\t" + \
		str(int(baltimore[b])) + " matches\n")

print("\n")

for br in broad_balt:
	print(str(br).split("\n")[0] + "\t" + \
	str(float(broad_balt[br]/num_line)*100).split("\n")[0] + "%" + "\t" + \
	str(int(broad_balt[br])) + " matches")
	if write:
		virus_broad.write(str(br).split("\n")[0] + "\t" + \
		str(float(broad_balt[br]/num_line)*100).split("\n")[0] + "%" + "\t" + \
		str(int(broad_balt[br])) + " matches\n")

if write:
	virus_family.close()
	virus_baltimore.close()
	virus_broad.close()

#print("\n" + "Pipe output to 'sort -k3nr' to rank from highest to lowest proportion")
