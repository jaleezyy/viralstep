#!/bin/bash

#set -e # exit if pipeline returns non-zero status
#set -u # unset variables == error
set -o pipefail # return value of last command to exit with non-zero status

# Goal is to generate a simple .csv database where each accession on NCBI Virus has its taxID and virus family
# Allows for built-in accessionID to taxID without relying on other tools
# Resulting output file can be shared to minimize requirements for running the baits pipeline
# Script will take all genbank and refseq accessions and use KronaTools (ktGetTaxIDfromAcc) to generate taxIDs
# Taxonkit will be used to generate lineage information (i.e. virus families) from the list of TaxIDs
# Merger scripts will combine the resulting output files of the above into <accessionID, taxID, virus family>

# TODO: Replace KronaTools requirement (hopefully faster than my previous attempts)
# TODO: VHost-Classifer or use database file provided myself for added host information (NCBI alternative?) - want to expand such that we can filter sequences by host (ex. mammalian viruses)
# TODO: Create quick parser to control for input

HELP="""
Wrapper script that aims to take a FASTA file and generate a classification database file (.csv) containing the following:

accessionID, taxonomyID, balitmore classifier (if reference file found), virus family

Arguments are as followed:
generate_db.sh -i <input_fasta> -o <output_dirname>

Output directory will contain intermediate files and will be found in $(pwd)
By default, the output directory will be called 'db_intermediate'

"""

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
input=0
output=0

while getopts ":i:o:" option; do
	case "${option}" in
		i) input=$OPTARG;;
		o) output=$OPTARG;;
	esac
done

if [ $input = 0 ]; then 
	echo "Missing input FASTA file (and path)"
	echo "$HELP"
	exit 1
fi

if [ $output = 0 ]; then
	echo "Output directory will be: db_intermediate"
	output="db_intermediate"
else
	echo "Output directory will be: ${output}"
fi

if [ -f ""$input"" ]; then
	#input_file=$(basename $input | rev | cut -d. -f1 | rev)
	input_file=${input%.*}
	current_dir=$(pwd)
else
	echo "Missing input FASTA file (and path)"
	echo "$HELP"
	exit 1
fi

# make directories
mkdir $current_dir/$output
mkdir $current_dir/$output/logs
num_seq=$(grep -c ">" $input)

# Download all known virus sequences or accessions
# ncbi-genome-download # genbank
# ncbi-genome-download # refseq
# TODO: add merge to generate 3 files: one genbank, one refseq, and one merged

checkpoint (){

# $1 = intermediate file to verify
# $2 = step number
	num_checkpoint=$(wc -l $1 | cut -d" " -f1)
	if [ $num_checkpoint -ne $num_seq ]; then
	echo -e "Something went wrong at (${2}/5) - check intermediate file: ${1}"
	clean_up
	fi

}

clean_up () {

	echo -e "\nCleaning up"
	mv accessionID_${input_file}.log $current_dir/$output/logs
	mv family_lineage_${input_file}.log $current_dir/$output/logs
	for file in *${input_file}*; do
		if [ ""$file"" != ""$input"" ]; then
			if [ ""$file"" != "db_${input_file}.csv" ] && [ ! -d ""$file"" ]; then
				mv $file $current_dir/$output
			fi
		fi
	done
	
	exit 0
}

# Isolate accession IDs (output will be where THIS script is run)
echo "(1/5) Extracting Accession IDs"
$script_dir/accessionID_v2.sh -s ${input}

checkpoint accessionID_${input_file}.txt 1

# Convert accession IDs to taxIDs and pipe through Taxonkit for lineage information
#cat accessionID_${input_file} | ktGetTaxIDfromAcc > taxID_${input_file}.txt
echo "(2/5) Converting to Taxonomy IDs and determining lineage"
cat accessionID_${input_file}.txt | ktGetTaxIDFromAcc | taxonkit lineage > lineage_${input_file}.taxonkit

checkpoint lineage_${input_file}.taxonkit 2

# Isolate virus family per taxID
echo "(3/5) Extracting virus family per sequence"
$script_dir/virus_family.sh -s lineage_${input_file}.taxonkit

checkpoint family_lineage_${input_file}.txt 3

# Merge above accessionIDs with taxIDs,virus_families
echo "(4/5) Creating database file"
ref=${script_dir}/reference_files/baltimore_ref.txt
$script_dir/merge_classify.py accessionID_${input_file}.txt family_lineage_${input_file}.txt db_${input_file}.csv $ref

checkpoint db_${input_file}.csv 4

# Calculate proportions
echo "(5/5) Calculating proportions of sequences per virus family"
# $script_dir/proportion.py family_lineage_${input_file}.txt | sort -k3nr > proportion_${input_file}.txt
$script_dir/proportion.py db_${input_file}.csv ${input_file}

# Move log files
echo -e "\nDatabase file created: db_${input_file}.csv"
clean_up
