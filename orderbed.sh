#!/bin/bash

## orderbed.sh


if [ $# -lt 2 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2017"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tOrder bed file by provided chromosome order.\n"
	echo -e "Usage:  \t`basename $0` [bed] [chrom] [--no_error] > [stdout]\n"
	echo -e "        \t[chrom]: A list of chromosomes. Only field 1 will be used. Use chromcheck.sh to get a chromosome list from a bed file."
	echo -e "        \t[--no_error]: If set, errors will not be produced even if only a part of chroms in bed file exist in the reference chromosome file.\n"
	exit 
fi

# 1) Files and variables

	bed=$1
	chrom=$2
	
	indf=$(head -n1 $bed | awk '{print NF+1}')
	cutf=$(head -n1 $bed | awk '{print NF}')
	
	if [ -n "$3" ] && [ $3 == "--no_error" ]
	then
		noerror=1
	else
		noerror=0
	fi

# 2) Processing
	
	## 2.1) Add chrom_index to chrom list
	cut -f 1 $chrom |
	awk '{OFS="\t"; print $1, NR}' |
	 
	## 2.2) Add defined chrom_index to bed file 
	awk -F "\t" '
	BEGIN{error_check='$noerror';}
	(NR==FNR){chrom_index[$1]=$2; next;} 
	{ if($1 in chrom_index){ print $0 "\t" chrom_index[$1]; next;}
		else if(error_check == 1) {next;}
		else { print "ERROR:" $1 " in bed not found in reference chrom list."; next;}
	}' - $bed |
	
	## 2.3) Sort by chrom_index, and then by start position 
	sort -k ''$indf','$indf'n' -k '2,2n' |
	cut -f '1-'$cutf''
	
exit;

	
###########################################
# Copyright (c) 2017, Kohta Ikegami
# All rights reserved.
# contact: ikgmk@uchicago.edu
###########################################	
