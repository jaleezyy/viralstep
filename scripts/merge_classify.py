#!/bin/env python

### Goal is to merge accessionID values with their corresponding taxID, balitmore classifier, virus_family

import sys
from collections import defaultdict

try:
	if "help" in sys.argv: raise IndexError
	input_ids = sys.argv[1]
	input_families = sys.argv[2]
	output = str(sys.argv[3])
	baltimore_ref = sys.argv[4]
except IndexError:
	print("""
Arguments are as followed:
	merge_classify.py <input_file_ids> <input_file_families> <output_file> <baltimore_ref_file>
	""")
	exit()

baltimore2family = defaultdict(list)
family2baltimore = {}

for line in open(baltimore_ref):
	family = line.split("\t")[0]
	classify = line.strip("\n").split("\t")[1] # baltimore classification (i.e. ssRNA, dsDNA, etc.)
	baltimore2family[classify].append(family)
	family2baltimore[family] = classify

out = open(output, "w+")

with open(input_ids) as f:
	with open(input_families) as f2:
		# line1 = AccID
		# line2 = taxID, virus_family
		# output = line1, line2[0], baltimore, line2[1] 
		for line1,line2 in zip(f,f2):
			accID = line1.strip("\n")
			id_and_family = line2.strip("\n")
			try:
				out.write(accID + "," + \
					id_and_family.split(",")[0].strip() + "," + \
					family2baltimore[id_and_family.split(",")[1]].strip() + "," + \
					id_and_family.split(",")[1].strip() + "\n")
			except (KeyError, IndexError): # write in "Unknown" for baltimore
				out.write(accID + "," + \
					id_and_family.split(",")[0].strip() + "," + \
					"Unknown," + \
					id_and_family.split(",")[1].strip() + "\n")

out.close()				
print("Done.")
exit()
			