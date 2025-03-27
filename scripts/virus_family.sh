#!/bin/sh

# From Taxonkit lineage output, isolate the TaxID, Virus family (ends in "\w*dae") for each line in the lineage data

# set -e # exit if pipeline returns non-zero status
# set -u # unset variables == error
# set -o pipefail # return value of last command to exit with non-zero status

# -s: SOURCE = lineage data (.txt)

while getopts ":s:" option; do
	case "${option}" in
		s) SOURCE=$OPTARG;;
		
esac
done

# Get filename
# count=$(echo $SOURCE | tr -cd '/' | wc -c)
# separation=$(($count+1))
# filename=$(echo $SOURCE | awk -F "/" '{print $var}' var="${separation}")
filename=${SOURCE%.*}

# Start of program
echo "Started at "$(date)"." > "family_"$filename".log" && echo "Started at "$(date)"."


# Check to make sure we only work with Virus lineage data (may be redundant)
# grep "Viruses" $SOURCE > "family_"$filename".txt_intermediate"

count=0
total=$(wc -l $SOURCE | awk -F " " '{print $1}')
echo $count/$total >> "family_"$filename".log" && echo $count/$total

while read line; do
	echo $(echo $line | awk -F " " '{print $1}')","$(echo $line | grep -o "\w*dae" || echo "Unclassified") | awk -F " " '{print $1}' >> "family_"$filename".txt"
	count=$(($count+1))
	echo $count/$total >> "family_"$filename".log" && echo $count/$total
	
done < $SOURCE


#rm "family_"$filename".txt_intermediate"
#awk -F " " '{print $1}' "family_"$SOURCE".txt_intermediate" > "family_"$SOURCE".txt"

# End of program
echo "Done at "$(date)"." >> "family_"$filename".log" && echo "Done at "$(date)"."

