#!/usr/bin/env bash

FASTA=0
THREADS=4

HELP="""
Executable will run through all the steps of ViralSTEP. Plan is to convert this script into the main ViralSTEP executable

Command:
	multi_viralstep.sh -f <fasta_file>

Flags:
	-f  :  Input FASTA file that will be used as reference material for ViralSTEP 
	-t  :  Number of threads. Default = 4
"""

while getopts "f:t:" option; do
	case "${option}" in
		f) FASTA=$OPTARG;;
		t) THREADS=$OPTARG;;
	esac
done

if [ $FASTA = 0 ] || [ ! -f $FASTA ];  then
	echo "Specify reference FASTA file for bait design"
	echo "$HELP"
	exit 1
fi

### Setup parameters
start=$PWD
fasta_name=$(basename -s .fa $(basename -s .fasta $FASTA))
fasta_path=$(dirname $FASTA)

### Create output directories for baits
out_dir=${fasta_path}/${fasta_name}_baits
mkdir -p $out_dir

# Tilebaits
viralstep tilebaits --default --outprefix ${fasta_name} --outdir ${out_dir} -i $FASTA

expected_output=${out_dir}/${fasta_name}-filtered-baits.fa

if [ ! -f $expected_output ]; then
	echo "Failed ViralSTEP"
	touch $expected_output
fi

# Melting temperature (default 65-95 degC)
tm_dir=${out_dir}/tm
mkdir -p $tm_dir
ln -s $(realpath $expected_output) $tm_dir/$(basename $expected_output)

cd $tm_dir
viralstep tmbaits -i $(basename $expected_output)
cd - # outside out_dir (i.e., ./out_dir/...)

expected_output=${tm_dir}/TmFiltered_${fasta_name}-filtered-baits.fa

if [ ! -f $expected_output ]; then
	echo "Failed melting temperature filter"
	touch $expected_output
fi

# Cross-hybridization
cross_dir=${out_dir}/cross
mkdir -p $cross_dir

viralstep crosshybaits -i $expected_output -o ${cross_dir}/crosshyb_tm_baits.fa -k Viruses --validate -t $THREADS

expected_output=${cross_dir}/crosshyb_tm_baits.fa

if [ ! -f $expected_output ]; then
	echo "Failed cross-hybridization filter"
	touch $expected_output
fi

# Redundancy (60% identity)
redund_dir=${out_dir}/redund
mkdir -p $redund_dir

viralstep redundbaits -i $expected_output -m 60 -o ${redund_dir}/redund_crosshyb_tm_baits.fa -t $THREADS

expected_output=${redund_dir}/redund_crosshyb_tm_baits.fa

if [ ! -f $expected_output ]; then
	echo "Failed redundancy filter"
	touch $expected_output
fi

num_baits=$(grep -c ">" $expected_output)

echo -e "Complete! ${num_baits} bait sequences proposed!"
exit 0
