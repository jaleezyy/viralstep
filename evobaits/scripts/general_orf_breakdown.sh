#!/usr/bin/env bash

input_dir=0
field="1"

### Parser
HELP="""
Usage:
bash general_orf_breakdown.sh -d <directory>

This script aims to produce a series of statistics per ORF for a given directory representing a single virus family. All output will be sent to stdout.

Flags:
	-d  :  Input directory representing a single virus family. Directory should be labelled as the virus family.
	-f  :  Field position where ORF name can be found across filenames (default = 1).

Credits: Script developed by Jalees A. Nasir, McArthur Lab, @Jaleezyy, 2023
Version: 0.2
"""

while getopts ":d:f:" option; do
	case "${option}" in
		d) input_dir=$OPTARG;;
		f) field=$OPTARG;;
	esac
done

if [[ $input_dir == 0 ]]; then
	echo "Missing input!"
	echo "$HELP"
	exit 1
elif [ ! -d $input_dir ]; then
	echo "Input must be a directory!"
	echo "$HELP"
	exit 1
fi

# store cut command to extract ORF name
cut_cmd="cut -d_ -f${field}"

# prioritize FASTA then DOT
# Required stats: # clusters, # clusters with 1 ORF
# separate analysis? virus families total, # virus families with 1 ORF
virus_family=$(basename $input_dir) 
echo -e "${virus_family}:"

clusters=0
single_clusters=0

# within each ORF, derive some stats
# number of sequences, average length of ORF sequence
echo -e "ORF\tNumber of Sequences\tNumber of Unique Virus Families"
for file in $(find $input_dir -name "*.fasta" ! -name "*unclustered.*" -print); do 
	# pull name
	orf=$(basename $file | cut -d. -f1 | $cut_cmd | awk '{print tolower($0)}')

	# Count number of sequences and number of clusters with only one sequence
	count=$(grep -c ">" $file)
	clusters=$(($clusters+1))
	if [[ $(grep -c ">" $file) -eq 1 ]]; then
		single_clusters=$(($single_clusters+1))
	fi
	
	# Count number of unique virus families
	# Headers expected to be {ref|id)_*viridae
	unique_vf=$(grep ">" $file | rev | $cut_cmd | rev | awk '{print tolower($0)}' | sort -u | wc -l)
	
	# report to stdout
	#echo -e "${orf}: ${count} sequence(s) (Number of unique virus families: ${unique_vf})"
	echo -e "${orf}\t${count}\t${unique_vf}"
done

echo -e ""
echo -e "Number of clusters identifed for ${virus_family}: ${clusters}"
echo -e "Number of clusters containing only one sequence within ${virus_family}: ${single_clusters}"
echo -e ""
echo -e "----------"
echo -e ""
