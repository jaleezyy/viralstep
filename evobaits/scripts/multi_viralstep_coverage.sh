#!/usr/bin/env bash

FASTA=0
THREADS=4
REF="/workspace/lab/mcarthurlab/nasirja/workspace/baits_pipeline/evo_baits/test_run/all/hyphy_single/spike_cluster.fasta" # set as default for now
DB="/workspace/lab/mcarthurlab/nasirja/workspace/baits_pipeline/evo_baits/test_run/all/hyphy_single/spike_cluster.fasta.csv"
TERM="Spike_cluster"

HELP="""
Executable will run through the coverage calculation portion of ViralSTEP. Plan is to convert this script into the main ViralSTEP executable

Command:
	multi_viralstep_coverage.sh -f <fasta_file> -r <fasta_ref> -d <fasta_ref_db>

Flags:
	-f  :  Input FASTA file that will be used as reference material for ViralSTEP. FASTA is a bait set, intermediate or final. 
	-r  :  Input FASTA file for bait target reference sequences. Represents all potential targets for a given bait set. Default is '${REF}'
	-d  :  Input FASTA ViralSTEP DB corresponding to bait target reference. Default is '${DB}'
        -s  :  Search term, be it virus family, Spike_cluster, or syndromic term. Default is '${TERM}'
	-t  :  Number of threads. Default = 4
"""

while getopts "f:t:r:d:s:" option; do
	case "${option}" in
		f) FASTA=$OPTARG;;
		r) REF=$OPTARG;;
		d) DB=$OPTARG;;
                s) TERM=$OPTARG;;
		t) THREADS=$OPTARG;;
	esac
done

if [ $FASTA = 0 ] || [ ! -f $FASTA ];  then
	echo "Specify reference FASTA file for bait design"
	echo "$HELP"
	exit 1
fi

if [ ! -f $REF ]; then
	echo -e "Missing FASTA reference ${REF}!"
	echo "$HELP"
	exit 1
fi

if [ ! -f $DB ]; then
	echo -e "Missing FASTA DB file ${DB}!"
	echo "$HELP"
	exit 1
fi

### Setup parameters
start=$PWD
fasta_file=$(basename $FASTA)
fasta_name=$(basename -s .fa $(basename -s .fasta $FASTA))
fasta_path=$(dirname $FASTA)
ref_path=$(realpath $REF)
db_path=$(realpath $DB)

cd $fasta_path

outputdir=${fasta_file}_coveragealign
outfile=${fasta_file}_coveragealign/merged_coverage.txt
outstats=${fasta_name}_plots
outstatsfile=${fasta_name}_plots/summary.txt

if [ ! -d $outputdir ] && [ ! -f $outfile ]; then
	viralstep coverage -i $fasta_file -r $ref_path -d $db_path
else
	#echo "Coverage files found! Skipping..."
	rm -r $outputdir
	viralstep coverage -i $fasta_file -r $ref_path -d $db_path
fi

if [ ! -d $outstats ] && [ ! -f $outstatsfile ]; then
	viralstep plots -o $outstats -c $outfile -s $TERM
else
	# echo "Stats files found! Skipping..."
	rm -r $outstats
	viralstep plots -o $outstats -c $outfile -s $TERM
fi

cov=$(grep "Average fold coverage for" ${outstatsfile} | cut -d " " -f 9)
if [[ $cov = '' ]]; then
	cov=0
fi
com=$(grep "Average % genome completeness for" ${outstatsfile} | cut -d" " -f 10)
if [[ $com = '' ]]; then
	com=0
fi 

echo -e "Coverage: ${cov}"
echo -e "Completeness: ${com}"

cd $start

exit 0

