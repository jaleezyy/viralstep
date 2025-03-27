#!/usr/bin/env bash

input_dir=0
skip_dot='false'
respir='false'
field="1"
warn=0
respir_families=('adenoviridae' 'orthomyxoviridae' 'anelloviridae' 'paramyxoviridae' 'arenaviridae' 'parvoviridae' 'bornaviridae' 'peribunyaviridae' 'coronaviridae' 'phenuiviridae' 'filoviridae' 'picornaviridae' 'reoviridae' 'flaviviridae' 'pneumoviridae' 'hantaviridae' 'polyomaviridae' 'herpesviridae' 'togaviridae' 'nairoviridae')
SCRIPT_PATH=$(dirname -- "${BASH_SOURCE[0]}")

### Parser
HELP="""
Usage:
bash general_cluster_breakdown.sh -d <directory>

This script aims to produce a series of statistics per virus family for a given directory representing a collection of virus family directories. A secondary script (generate_orf_breakdown.sh) will be invoked to analyse each individual virus family directory after summary statistics across all.

Flags:
	-d  :  Input directory representing a collection of virus family clusters.
	-f  :  Field position where ORF name can be found across filenames (default = 1).
	-s  :  Skip searching for DOT output. Rely on FASTA input only.
	-r  :  Filter results for 'respiratory' TODO: make this a search term.

Credits: Script developed by Jalees A. Nasir, McArthur Lab, @Jaleezyy, 2023
Version: 0.2
"""

while getopts ":d:f:sr" option; do
	case "${option}" in
		d) input_dir=$OPTARG;;
		f) field=$OPTARG;;
		s) skip_dot='true';;
		r) respir='true';;
	esac
done

if [ -f $SCRIPT_PATH/general_orf_breakdown.sh ]; then
	cscript=$SCRIPT_PATH/general_orf_breakdown.sh
else
	cscript=0
fi

echo $cscript

if [[ $input_dir == 0 ]]; then
	echo "Missing input!"
	echo "$HELP"
	exit 1
elif [ ! -d $input_dir ]; then
	echo "Input must be a directory!"
	echo "$HELP"
	exit 1
fi

# manual count of clusters slower, but can allow potential discrepencies to be detected per ORF cluster
total_clusters=0
vf=0
single_vf=0

# count the number of directories
# search within each directory to count the number of ORF clusters, highlight single ones

#count=$(ls $input_dir | wc -l)
for dir in $input_dir/*; do
	if [ -d $dir ]; then
		vf_name=$(basename $dir | cut -d. -f1 | cut -d_ -f1 | awk '{print tolower($0)}')
		# check if we're only filtering for respir, continue if vf not a match
		if [ $respir == 'true' ] && [[ ! " ${respir_families[*]} " =~ [[:space:]]${vf_name}[[:space:]] ]]; then
			continue
		fi
		
		vf=$(($vf+1))
		
		# when counting total number of clusters, check if we're only using FASTA input
		if [[ ! $cscript == 0 ]] && ([ $(find $dir -name "*.fasta" ! -name "*unclustered.*" -print | wc -l) -gt 0 ] || [ $(find $dir "*.dot" ! -name "*unclustered.*" -print | wc -l) -gt 0 ]); then
			count_fasta=$(find $dir -name "*.fasta" ! -name "*unclustered.*" -print | wc -l)
			count_dot=$(find $dir -name "*.dot" ! -name "*unclustered.*" -print | wc -l)
			
			if [ $skip_dot != 'true' ]; then
				# if number of ORFs is 1, flag it
				if ([ $count_fasta -eq 1 ] || [ $count_dot -eq 1 ]); then
					single_vf=$(($single_vf+1))
					total_clusters=$(($total_clusters+1))
				elif [ $count_fasta -eq $count_dot ]; then
					total_clusters=$(($total_clusters+$count_fasta)) # fasta == dot
				else
					if [ $count_fasta -lt $count_dot ]; then
						total_clusters=$(($total_clusters+$count_fasta)) # dot > fasta
						warn=1
					else
						total_clusters=$(($total_clusters+$count_dot)) # fasta > dot
						warn=1
					fi
				fi
			# using only FASTA output
			else
				if [ $count_fasta -eq 1 ]; then
					single_vf=$(($single_vf+1))
				fi
				total_clusters=$(($total_clusters+$count_fasta))
			fi
			bash $cscript -d $dir -f $field
		fi
	fi
done

echo -e ""
echo -e "Total number of clusters found: ${total_clusters}"
echo -e "Total number of virus families: ${vf}"
echo -e "Total number of virus families with only one ORF: ${single_vf}"
echo -e ""

if [[ $warn != 0 ]]; then
	echo -e "WARNING: Some clusters could not be accurately counted! The lowest values were counted so there may be a possible discrepency in the total number of clusters!"
	echo -e "Re-running with '-s' to skip DOT file assessment in cluster counts and use FASTA files only will help with consistency!"
fi
