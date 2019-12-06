#!/bin/bash

## chromcheck.sh

if [ $# -lt 2 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2017"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tPrint the order of chromosomes in a file.\n"
	echo -e "Usage:  \t`basename $0` [bed] [inspect every n lines] > [stdout]"
	echo -e "\n"
	exit 
fi

infile=$1
res=$2

sed -n '0~'$res'p' $infile | cut -f1 | uniq

exit

###########################################
# Copyright (c) 2017, Kohta Ikegami
# All rights reserved.
# contact: ikgmk@uchicago.edu
###########################################	