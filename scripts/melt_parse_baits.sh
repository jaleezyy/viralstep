#!/bin/env bash

# Run melt.pl and parse the data such that only Tm data outputted to a new file
# melt.pl run with given probe .fasta
# parsing occurs such that the data is first split between the hybrid software first run in melt.pl and the actual d* and Tm data
# The 4th column contains data for Tm which is separated and the header is removed
# An output file named with Tm_$SOURCE.txt is generated in the current directory

#set -e # exit if pipeline returns non-zero status
set -u # unset variables == error
set -o pipefail # return value of last command to exit with non-zero status

# -s: SOURCE = .fasta file with probe sequences 

while getopts ":s:o:" option; do
	case "${option}" in
		s) SOURCE=$OPTARG;;
		o) OUTPUT=$OPTARG;;
esac
done

source_name=$(basename $SOURCE)".symlink"

# generate symlink to SOURCE (will allow script to run from anywhere)
echo "(1/3) Linking input FASTA"
ln -s ${SOURCE} ${source_name}

echo "(2/3) Calculating melting temperatures"
melt.pl -n RNA -t 65 -C 1.89e-9 ${source_name} > "${OUTPUT}_intermediate"

count=$(grep -n "dG" ${OUTPUT}"_intermediate" | awk -F ":" '{print $1}') # find line number where data splits between both halves of melt.pl
total=$(wc -l ${OUTPUT}"_intermediate" | awk -F " " '{print $1}') # get total number of lines so a range can be set

# parse data to isolate Tm data
echo "(3/3) Parsing melting temperature"
echo "Melting Temperature (C)" > $OUTPUT
cat "${OUTPUT}_intermediate" | sed -n $count,$total'p' | awk '{print $4}' | tail -$(($count-1)) | head -$(($count-1)) >> $OUTPUT

# melt.pl specific output (removal)
#rm ${source_name}{".65.ext",".65.plot",".ct",".dG",".run"}
# echo "(4/4) Cleaning up"
rm "${OUTPUT}_intermediate"
rm ${source_name}
