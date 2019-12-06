#!/bin/bash

## mergebg.sh


if [ $# -lt 1 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2017"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tMerge neighboring bed regions with identical scores.\n"
	echo -e "Usage:  \t`basename $0` [bed] > [stdout]\n"

	exit 
fi

# 1) Files and variables
	bed=$1
	
# 2) Processing

	awk 'BEGIN {chr=0} 
		{	if(chr == 0) {chr=$1; start=$2; stop=$3; score=$4; next;}
			else if ($1== chr && $2 == stop && $4 == score) {stop=$3; next;}
			else {OFS="\t"; print chr, start, stop, score; chr=$1; start=$2; stop=$3; score=$4; next;}
		}
		END {OFS="\t"; print chr, start, stop, score;}' $bed
exit;
