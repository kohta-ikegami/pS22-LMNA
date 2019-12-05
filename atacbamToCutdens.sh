#!/bin/bash

## atacbamToCutdens.sh
## Requirement: bedtools v2.25.0

if [ $# -lt 6 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2016"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tGenerate Tn5 cut density profile.\n"
	echo -e "Usage:  \t`basename $0` [.bam (sorted)] [chrom] [read LEN] [extension] [scale] [out seed]"
	echo -e "\n"
	exit 
fi

# Files and parameters
	bam=$1
	genome=$2
	readLen=$3
	ext=$4
	scale=$5
	seed=$6
	bgout="$seed".bg
	bwout="$seed".bw
	cigar="$readLen"M	
		
	# Cut position and extension
	## Forward strand: cut center, bamCol4-1+4;
	## Reverse strand: cut center, bamCol4-1+readLen-5

	fStart=$((-1+4-$ext))				# Number to be added to bamCol4.
	fEnd=$((-1+4+$ext))					# Number to be added to bamCol4.
	rStart=$((-1+$readLen-5-$ext))		# Number to be added to bamCol4.
	rEnd=$((-1+$readLen-5+$ext))		# Number to be added to bamCol4.

# Processing
	# 1) Make bedGraph
	# Forward strand FLAGs
		# 99:	read paired, read mapped in proper pair, mate reverse strand, first in pair
		# 163:	read paired, read mapped in proper pair, mate reverse strand, second in pair
	# Reverse strnd FLAGs
		# 83:	read paired, read mapped in proper pair, read reverse strand, first in pair
		# 147:	read paired, read mapped in proper pair, read reverse strand, second in pair

	samtools view $bam |
	awk '($6=="'$cigar'") { 
		if ($2==99 || $2==163) {chr=$3; start=$4+'$fStart'; stop=$4+'$fEnd'}  
		else if ($2==83 || $2==147) {chr=$3; start=$4+'$rStart'; stop=$4+'$rEnd'}
		if (start<0){start=0;}
		OFS="\t"; print chr, start, stop; next;}' |
	bedtools genomecov -i stdin -bg -g $genome -scale $scale > $bgout

	# 2) Make bigwig
	bedGraphToBigWig $bgout $genome $bwout
	
exit
	
###########################################
# Copyright (c) 2016, Kohta Ikegami
# All rights reserved.
# contact: ikgmk@uchicago.edu
###########################################	