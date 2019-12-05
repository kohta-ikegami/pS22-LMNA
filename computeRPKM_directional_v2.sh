#!/bin/bash

## computeRPKM_directional.sh

if [ $# -lt 7 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2017"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tCompute RPKM from directional RNA-seq data.\n"
	echo -e "Usage:  \t`basename $0` [bam] [seed] [-sorted/-unsorted] [transcript] [cpu] [memory] [rpkm dir]\n"

	echo -e "        \t[bam]              \taccepted_hits.(sorted.)bam from tophat."
	echo -e "        \t[seed]             \tseed.rpkm and seed.count will be generated."	
	echo -e "        \t[-sorted/-unsorted]\tIs Bam file sorted?"
	echo -e "        \t[transcript]       \t*PerGeneSymbolTranscript.bed made by makePerGeneSymbolTranscriptBed.sh"
	echo -e "        \t[cpu]              \tNumber of threads to be used for samtools sort. Only with -unsorted option."
	echo -e "        \t[memory]           \tMemory per thread. e.g. 5 (Request in pbs [cpu] x [mem])."
	echo -e "        \t[rpkm dir]         \tDirectory for seed.rpkm.\n"
	
	echo -e "        \tTemp files will be generated in current working directory.\n"
		
	exit 
fi

# History
# v1.0:
# v2.0: Use scratch space for temp files 

# 1) Variables

	bam=$1		# tophat bam file
	seed=$2
	sortKey=$3	# bamfile is sorted or not 
	transcript=$4		#  chrom, txStart, txEnd, transcript.id, length, strand, gene.symbol, geneType, exonCount, exonStarts, exonEnds 
	mapq=50		# MAPQ
	cpu=$5
	mem=$6
	rpkmdir=$7
	rpkm="$seed".rpkm
	
# 2) Temp file names
	
	myram=$seed.$RANDOM.$RANDOM.$RANDOM
	
	chromtemp=/scratch/kikegami/chromtemp.$myram.temp
	
	plusexon=/scratch/kikegami/plusexon.$myram.temp
	minusexon=/scratch/kikegami/minusexon.$myram.temp
	
	plusbg=/scratch/kikegami/plusbg.$myram.temp
	minusbg=/scratch/kikegami/minusbg.$myram.temp

	plusexoncov=/scratch/kikegami/plusexoncov.$myram.temp
	minusexoncov=/scratch/kikegami/minusexoncov.$myram.temp
		
		
# 3) Make temp genome file from bam
	module load gcc/6.2.0 bedtools/2.25.0 samtools/1.3.1
	
	bamToChromfile.sh $bam > $chromtemp
	echo -e "Made chromtemp file."

# 4) Make exon file for each strand
	# Remove header
	sed '1d' $transcript | 		
	# Expand columns 10 (exon starts) and 11 (exon ends)	
	bedtools expand -c 10,11 |		
	# For each strand, report chr, exonStart, exonEnd, transcript.id (NOT SORTED)
	awk '{ if($6 == "+") {OFS="\t"; print $1, $10, $11, $4 > "'$plusexon'"} else {OFS="\t"; print $1, $10, $11, $4 > "'$minusexon'"} }'

# 5) Set variables for RPKM calculation
	# Comma separated field numbers to group (all fileds in the target bed) (for STEP2)	
	grp=$(head -n 1 $plusexon | awk '{ for (i=1; i<NF; i++) {printf "%s,", i;} } END {printf i}')
	echo -e "-grp columns:" $grp
	
	# The number of the field for which sum is computed (the "-c" column for groupby) (for STEP2)
	# NF + bedgraph (4) + intersectSize (1) + intersectSumPerbasecov (1) => NF+6
	cfield=$(head -n 1 $plusexon | awk '{print NF+6}')
	echo -e "-c columns": $cfield
	
	# Read length (for STEP3)
	rlength=$(samtools view -q $mapq $bam | head -n 1 | awk '{print length($10)}')	
	echo -e "Read length:" $rlength
	
	# Total number of reads (for STEP3)
	readnum=$(samtools view -c -q $mapq $bam)
	echo -e "Total number of reads in "$bam" at MAPQ "$mapq":" $readnum

# 6) Make bedgraph for each strand
	module load gcc/6.2.0 bedtools/2.26.0 samtools/1.3.1
	
	if [ $sortKey = "-sorted" ]; then
		samtools view -u -q $mapq $bam | 
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand + > $plusbg
		
		samtools view -u -q $mapq $bam | 
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand - > $minusbg
		
    elif [ $sortKey = "-unsorted" ]; then
    	samtools view -uhq $mapq $bam | 
		samtools sort -@ $cpu -m $mem -T $RANDOM.$RANDOM |  
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand + > $plusbg
		
    	samtools view -uhq $mapq $bam | 
		samtools sort -@ $cpu -m $mem -T $RANDOM.$RANDOM |  
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand - > $minusbg
			    
	else 
		echo -e "\tSort option [-sorted/-unsorted] required. Exit."
		exit
	fi 
		# Output is [1-chr, 2-start, 3-stop, 4-basecoverage]

# 7) Compute exon base coverage for each strand
	module load gcc/6.2.0 bedtools/2.25.0 samtools/1.3.1
		
	# 4.1) Plus-strand transcripts (use plus-strand exon and minus-strand bg)
	## Compute sum of perbase coverage in each exon (-wao reports the original A and B entries plus the number of basepairs of overlap )
	bedtools sort -faidx $chromtemp -i $plusexon |	
	bedtools intersect -a stdin -b $minusbg -sorted -wao -g $chromtemp -nonamecheck |
	awk '{OFS="\t"; print $0, $(NF-1)*$(NF)}' |
	sort -f -k1,1 -k2,2n -k4,4 |
	bedtools groupby -i stdin -g $grp -c $cfield -o sum > $plusexoncov
		# Output is [1-chr, 2-exonStart, 3-exonStop, 4-transcript.id, 5-exonSumPerbasecov]
	
	# 4.2) Minus-strand transcripts (use minus-strand exon and plus-strand bg)
	## Compute sum of perbase coverage in each exon (-wao reports the original A and B entries plus the number of basepairs of overlap )
	bedtools sort -faidx $chromtemp -i $minusexon |	
	bedtools intersect -a stdin -b $plusbg -sorted -wao -g $chromtemp -nonamecheck |
	awk '{OFS="\t"; print $0, $(NF-1)*$(NF)}' |
	sort -f -k1,1 -k2,2n -k4,4 |
	bedtools groupby -i stdin -g $grp -c $cfield -o sum > $minusexoncov
		# Output is [1-chr, 2-exonStart, 3-exonStop, 4-transcript.id, 5-exonSumPerbasecov]	
	
# 8) Compute RPKM
	
	## Compute exon size for each exon
	cat $plusexoncov $minusexoncov |
	awk '{OFS="\t"; print $0, $3-$2}' |
		# Output is [1-chr, 2-exonStart, 3-exonStop, 4-transcript.id, 5-exonSumPerbasecov, 6-exonSize]
	
	## Sort exons by transcript.id
	sort -sf -k4,4 |
		# -s for stable sort. -f for ignoring case.
		 
	## For each transcript, get sum of the exon sum of perbasecoverage and exon size
	bedtools groupby -g 4 -c 5,6 -o sum,sum -full |
		# Output is [1-chr, 2-exonStart, 3-exonStop, 4-transcript.id, 5-exonSumPerbasecov, 6-exonSize, 7-transcriptSumPerbasecov, 8-transcriptSumExonSize]
	
	cut -f4,7,8 |
		# Output is [1-transcript.id, 2-transcriptSumPerbasecov, 3-transcriptSumExonSize]

	## For each transcript.id, compute RPKM (Reads Per Kilobase of transcript per Million mapped reads)
	### rpkm = [transcriptSumPerbasecov ($2)]/[read length ($readlength)]/[transcriptSumExonSize ($3/1000)]/[total num of aligned reads in million ($readnum/10^6))] 
	awk '{OFS="\t"; print $1, $3, $2, $2/'$rlength'/($3/1000)/('$readnum'/10^6)}' |
	
	## Sort by transcript.id and report
	sort -sf -k1,1 > $rpkmdir/$rpkm
	# Output is [1-transcript.id, 2-transcriptSumExonSize, 3-transcriptSumPerbasecov, 4-transcriptRPKM]	

# 9) Remove temp files
	rm $chromtemp $plusexon $minusexon $plusbg $minusbg $plusexoncov $minusexoncov


exit




###########################################
# Copyright (c) 2016, Kohta Ikegami
# All rights reserved.
# contact: ikgmk@uchicago.edu
###########################################	

