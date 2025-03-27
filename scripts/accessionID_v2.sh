#!/bin/sh


# Isolate headers from a FASTA file, and corresponding accessionID

#set -e # exit if pipeline returns non-zero status
#set -u # unset variables == error
#set -o pipefail # return value of last command to exit with non-zero status

# -s: SOURCE = .fasta sequence file

while getopts ":s:" option; do
	case "${option}" in
		s) SOURCE=$OPTARG;; # FASTA file with sequence header
		
esac
done

# Get filename
# count=$(echo $SOURCE | tr -cd '/' | wc -c)
# separation=$(($count+1))
# filename=$(echo $SOURCE | awk -F "/" '{print $var}' var="${separation}")
filename=${SOURCE%.*}

# Start of program
echo "Started at "$(date)"." > "accessionID_"$filename".log" && echo "Started at "$(date)"."

grep ">" $SOURCE > $filename".txt_intermediate" # isolate to create a list of headers with 1 per line

intermediate=$filename".txt_intermediate" # assign newly generated header file as intermediate

count=0
total=$(wc -l $intermediate | awk -F " " '{print $1}')
echo $count/$total >> "accessionID_"$filename".log" && echo $count/$total

while read line; do
	#echo $line | awk -F " " '{print $1}' | awk -F ">" '{print $2}' | awk -F "." '{print $1}' >> "accessionID_"$filename".txt"
	accID=$(echo $line | awk -F " " '{print $1}' | awk -F ">" '{print $2}') 
	if [ $(echo "${accID}" | tr -cd "_" | wc -c) -gt 0 ]; then
		echo "${accID}" | cut -d_ -f1,2 >> "accessionID_"$filename".txt"
	else
		echo "${accID}" | cut -d_ -f1 >> "accessionID_"$filename".txt" # accounts for bait headers
	fi
	count=$(($count+1))
	echo $count/$total >> "accessionID_"$filename".log" && echo $count/$total
	#echo $line >> "test.txt"
done < $intermediate

rm $intermediate

# End of program
echo "Done at "$(date)"." >> "accessionID_"$filename".log" && echo "Done at "$(date)"."

