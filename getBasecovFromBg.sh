#!/bin/bash

# getBasecovFromBg.sh

# Copyright (c) 2016, Kohta Ikegami
# All rights reserved.
# October 22, 2016
# 
if [ $# -lt 3 ] # '-lt' stands for less than 
then
	echo -e "\tby Kohta Ikegami 2015; all rights reserved"
  	echo -e "\n\tUsage: `basename $0` [target] [bedgraph] [chrom] | stdout"
  	echo -e "\n\t[target]\tBed file. For each interval, the sum of perbase count (or \"area\") will be computed."
  	echo -e "\t[bedgraph]\tBedgraph file of perbase count (not normalized)."
	echo -e "\t[chrom]\tChromosome length file. Typically, BAMsort version. Sort target by -k1,1V -k2,2n"
	exit 
fi

# 1) Modules
	# bedtools map in bedtools/2.26.0 has a bug. Therefore, this will use 2.25.0 
	module load bedtools/2.25.0
	
# 2) Variables
	afile=$1
	bedgraph=$2
	chrom=$3

# 3) Parameters in groupby
	## Comma separated field to group (all fileds in the target bed)	
	grp=$(head -n 1 $afile | awk '{ for (i=1; i<NF; i++) {printf "%s,", i;} } END {printf i}')

	## The number of the field for which sum is computed (the "-c" column for groupby)
	cfield=$(head -n 1 $afile | awk '{print NF+6}')

# 4) Compute sum of perbase coverage 
	bedtools intersect -a $afile -b $bedgraph -sorted -wao -g $chrom -nonamecheck | \
	awk '{OFS="\t"; print $0, $(NF-1)*$(NF)}'| \
	bedtools groupby -i stdin -g $grp -c $cfield -o sum

exit

