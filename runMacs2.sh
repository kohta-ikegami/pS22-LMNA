#!/bin/bash

## runMacs2.sh
## Requirement: macs2, bedtools v2.25.0, bedClip, bedGraphToBigWig

if [ $# -lt 7 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2017"
	echo -e "Version:\tv1.2"
	echo -e "Summary:\tCall peak and generate input-normalized signal profile in bg and bw.\n"
	echo -e "Usage:  \t`basename $0` [chip.bam] [input.bam] [species (hs/mm/ce)] [chrom] [ext] [dir] [seed]\n"
	echo -e "        \t- If multiple bam files, comma separate."
	echo -e "        \t- For seed, Preferred naming convention is [Ab_Cell_macs2hg19_tXXXtXXXcXXXcXXX]."
	echo -e "        \t- Request 4 cores, 10gb memory, and 8-hr walltime."
	echo -e "\n"
	exit 
fi

# Version changes
	# v.1.1
	# Added modules load In # Modules
	# Added sort -k1,1 --stable in # Processing Step 2.
	# Added "_k1sort" in # Files and parameters 
	# v1.2
	# Changed to be able to run on gardner
	
# Modules
#	module load MACS2/2.1.0 bedtools/2.25.0
	module load gcc/6.2.0
	module load bedtools/2.26.0 python/2.7.13 
	
# Files and parameters
	tbam=${1//,/ }
	cbam=${2//,/ }
	
	species=$3
	chrom=$4
	ext=$5
	dir=$6			
	seed=$7

	tbg="$seed"_treat_pileup.bdg
	cbg="$seed"_control_lambda.bdg
	unclipfebg="$seed"_fe.unclip

	febg="$seed"_fe_k1sort.bg
	febw="$seed"_fe_k1sort.bw

# Processing

# 1) MACS2
	echo -e `date` "\tStarting macs2 callpeak"
	macs2 callpeak -t $tbam -c $cbam -f BAM -g $species -n $seed --nomodel --extsize $ext --outdir $dir --SPMR --bdg --call-summits

# 2) Fold-enrichment bedgraph
	echo -e `date` "\tStarting bedtools unionbedg"
	bedtools unionbedg -i $tbg $cbg | 
	awk '{printf "%s\t%s\t%s\t%.3f\n", $1,$2,$3,$4/$5}' |
	sort -k1,1 --stable > $unclipfebg

# 3) Clip fold-enrichment bedgraph
	echo -e `date` "\tStarting bedClip"	
	bedClip $unclipfebg $chrom $febg

# 4) Convert to bigwig
	echo -e `date` "\tStarting bedGraphToBigWig"	
	bedGraphToBigWig $febg $chrom $febw

# 5) Gzip bg file
	gzip $febg

# 6) Remove temp files
	rm $tbg $cbg $unclipfebg

# 7) Software versions
	echo -e "Software versions: "
	macs2 --version
	bedtools --version

exit
