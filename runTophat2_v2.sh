#!/bin/bash

# runTophat2_v2.sh

if [ $# -lt 4 ]; then # '-lt' stands for less than

	echo -e "\n\tUsage: `basename $0` [fastq] [outdir] [cpu] [index seed]\n"
	
	echo -e "\t[fastq]         If more than one, separate by commas."	
	echo -e "\t[outdir]        Output directory in full path."
	echo -e "\t[cpu]           Num of CPU cores requested in pbs."
	echo -e "\t[index seed]    Full path with seed for bowtie index. e.g. /path/hg19ercc92\n"

	exit 
fi

# 1) Variables
	fastq=$1
	outdir=$2
	cpu=$3	
	index=$4
	
# 2) Outdir
	if [ ! -d $outdir ]; then mkdir $outdir; fi
	
# 3) Run tophat
	echo -e "Starting tophat2"
	
	module load gcc/6.2.0 tophat/2.1.1 bowtie2/2.1.0
	tophat -o $outdir -p $cpu --keep-fasta-order $index $fastq 
	
	echo -e "Finished running tophat."

exit

