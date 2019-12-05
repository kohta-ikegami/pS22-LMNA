#!/bin/bash

## Bruseq_pipeline_v1.sh

if [ $# -lt 8 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2017"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tCompute RPKM from directional Bru-seq data. The entire gene body will be used for RPKM calculation.\n"
	echo -e "Usage:  \t`basename $0` [bam] [seed] [-sorted/-unsorted] [transcript] [cpu] [memory] [rpkm dir] [-bgbw]\n"

	echo -e "        \t[bam]              \tBowtie-aligned bam. Reads with MAPQ >20 (-q 20) AND read-mapped-in-proper-pair (-f 3) will be used."
	echo -e "        \t[seed]             \tseed.rpkm and seed.count will be generated. Will also be used in bg/bw file names."	
	echo -e "        \t[-sorted/-unsorted]\tIs Bam file sorted?"
	echo -e "        \t[transcript]       \t*PerGeneSymbolTranscript.bed made by makePerGeneSymbolTranscriptBed.sh"
	echo -e "        \t[cpu]              \tNumber of threads to be used for samtools sort. Only with -unsorted option."
	echo -e "        \t[memory]           \tMemory per thread. e.g. 10 (Request in pbs [cpu] x [mem])."
	echo -e "        \t[rpkm dir]         \tDirectory for seed.rpkm."
	echo -e "        \t[-bgbw]            \tGenerate bedgraph and bigwig files.\n"
	
	echo -e "        \tTemp files will be generated in current working directory.\n"
		
	exit 
fi

# History
# v1.0: 
# v1.1: cleaned up for clarity. Exon name changed to txpt

# 1) Variables

	bam=$1
	seed=$2
	sortKey=$3	# bamfile is sorted or not 
	transcript=$4		#  chrom, txStart, txEnd, transcript.id, length, strand, gene.symbol, geneType, exonCount, exonStarts, exonEnds 
	mapq=20		# MAPQ
	flag=3
	cpu=$5
	mem=$6
	rpkmdir=$7
	rpkm="$seed".rpkm
	bgbw=$8
	
# 2) Temp file names
	chromtemp=chromtemp.$RANDOM.$RANDOM.temp
	
	plustxpt=plustxpt.$RANDOM.$RANDOM.temp
	minustxpt=minustxpt.$RANDOM.$RANDOM.temp
	
	plusbg=plusbg.$RANDOM.$RANDOM.temp
	minusbg=minusbg.$RANDOM.$RANDOM.temp

	plustxptcov=plustxptcov.$RANDOM.$RANDOM.temp
	minustxptcov=minustxptcov.$RANDOM.$RANDOM.temp
		
		
# 3) Make temp genome file from bam
	module load gcc/6.2.0 bedtools/2.25.0 samtools/1.3.1
	
	bamToChromfile.sh $bam > ./$chromtemp
	echo -e "Made chromtemp file."

# 4) Make exon file for each strand
	# In Bru-seq, do not "expand" exons. Instead, use the whole transcript as one exon.

	# Remove header
	sed '1d' $transcript | 		
	# For each strand, report chr, transcriptStart, transcriptEnd, transcript.id (NOT SORTED)
	awk '{ if($6 == "+") {OFS="\t"; print $1, $2, $3, $4 > "./'$plustxpt'"} else {OFS="\t"; print $1, $2, $3, $4 > "./'$minustxpt'"} }'

# 5) Set variables for RPKM calculation
	# Comma separated field numbers to group (all fileds in the target bed) (for STEP2)	
	grp=$(head -n 1 $plustxpt | awk '{ for (i=1; i<NF; i++) {printf "%s,", i;} } END {printf i}')
	echo -e "-grp columns:" $grp
	
	# The number of the field for which sum is computed (the "-c" column for groupby) (for STEP2)
	# NF + bedgraph (4) + intersectSize (1) + intersectSumPerbasecov (1) => NF+6
	cfield=$(head -n 1 $plustxpt | awk '{print NF+6}')
	echo -e "-c columns": $cfield
	
	# Total number of mapped bases (for STEP3) (=Sum of mapped fragment lengths. "-f 67" filters for first in pair AND read mapped in proper pair AND read paired. i.e. half of the reads.)
	totbase=$(samtools view -f 67 -q $mapq $bam | awk '{ if($2==99){sum += $9} else if ($2==83) {sum -= $9} else {print "error"} } END {print sum}')
	echo -e "Total number of mapped bases in "$bam" at MAPQ "$mapq" with FLAG "$flag":" $totbase

# 6) Make bedgraph for each strand
	module load gcc/6.2.0 bedtools/2.26.0 samtools/1.3.1
	
	if [ $sortKey = "-sorted" ]; then
		samtools view -u -q $mapq -f $flag $bam | 
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand + -pc > ./$plusbg
		
		samtools view -u -q $mapq -f $flag $bam | 
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand - -pc > ./$minusbg
		
    elif [ $sortKey = "-unsorted" ]; then
    	samtools view -uhq $mapq -f $flag $bam | 
		samtools sort -@ $cpu -m $mem -T $RANDOM.$RANDOM |  
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand + -pc > ./$plusbg
		
    	samtools view -uhq $mapq -f $flag $bam | 
		samtools sort -@ $cpu -m $mem -T $RANDOM.$RANDOM |  
		bedtools genomecov -ibam stdin -bg -g $chromtemp -split -strand - -pc > ./$minusbg
			    
	else 
		echo -e "\tSort option [-sorted/-unsorted] required. Exit."
		exit
	fi 
		# Output is [1-chr, 2-start, 3-stop, 4-basecoverage]

# 7) Compute transcript base coverage for each strand
	module load gcc/6.2.0 bedtools/2.25.0 samtools/1.3.1
		
	# 4.1) Plus-strand transcripts (use plus-strand transcript and minus-strand bg)
	## Compute sum of perbase coverage for each transcript (-wao reports the original A and B entries plus the number of basepairs of overlap )
	bedtools sort -faidx $chromtemp -i ./$plustxpt |	
	bedtools intersect -a stdin -b ./$minusbg -sorted -wao -g $chromtemp -nonamecheck |
	awk '{OFS="\t"; print $0, $(NF-1)*$(NF)}' |
	sort -f -k1,1 -k2,2n -k4,4 |
	bedtools groupby -i stdin -g $grp -c $cfield -o sum > ./$plustxptcov
		# Output is [1-chr, 2-transcriptStart, 3-transcriptStop, 4-transcript.id, 5-transcriptSumPerbasecov]
	
	# 4.2) Minus-strand transcripts (use minus-strand transcript and plus-strand bg)
	## Compute sum of perbase coverage in each transcript (-wao reports the original A and B entries plus the number of basepairs of overlap )
	bedtools sort -faidx $chromtemp -i ./$minustxpt |	
	bedtools intersect -a stdin -b ./$plusbg -sorted -wao -g $chromtemp -nonamecheck |
	awk '{OFS="\t"; print $0, $(NF-1)*$(NF)}' |
	sort -f -k1,1 -k2,2n -k4,4 |
	bedtools groupby -i stdin -g $grp -c $cfield -o sum > ./$minustxptcov
		# Output is [1-chr, 2-transcriptStart, 3-transcriptStop, 4-transcript.id, 5-transcriptSumPerbasecov]	
	
# 8) Compute RPKM
	
	## Compute transcript size
	cat $plustxptcov $minustxptcov |
	awk '{OFS="\t"; print $0, $3-$2}' |
		# Output is [1-chr, 2-transcriptStart, 3-transcriptStop, 4-transcript.id, 5-transcriptSumPerbasecov, 6-transcriptSize]
	
	## Sort transcripts by transcript.id
	sort -sf -k4,4 |
		# -s for stable sort. -f for ignoring case.
		 
	## For each transcript, get sum of the transcript sum of perbasecoverage (should be only one transcript per transcript.id) and transcript size
	bedtools groupby -g 4 -c 5,6 -o sum,sum -full |
		# Output is [1-chr, 2-transcriptStart, 3-transcriptStop, 4-transcript.id, 5-transcriptSumPerbasecov, 6-transcriptSize, 7-transcriptSumPerbasecov, 8-transcriptSumTranscriptSize]
	
	cut -f4,7,8 |
		# Output is [1-transcript.id, 2-transcriptSumPerbasecov, 3-transcriptSumTranscriptSize]

	## For each transcript.id, compute depth- and size-normalized coverage (Basecov Per Kilobase of transcript per total mapped bases)
	### rpkm = [transcriptSumPerbasecov ($2)]/[transcriptSumTranscriptSize ($3/1000)]/[total mapped bases ($totbase)] 
	awk '{OFS="\t"; print $1, $3, $2, $2/($3/1000)/'$totbase'}' |
	
	## Sort by transcript.id and report
	sort -sf -k1,1 > $rpkmdir/$rpkm
	# Output is [1-transcript.id, 2-transcriptSumTranscriptSize, 3-transcriptSumPerbasecov, 4-transcriptRPKM]	
	
	echo -e "Finished writing rpkm file."
	
# 9) Compute strand-specific normalized tracks

if [ $bgbw == "-bgbw" ]; then 
	# file names
	norm_plusbg="$seed"_np.bg
	norm_minusbg="$seed"_nm.bg
	norm_plusbw="$seed"_np.bw
	norm_minusbw="$seed"_nm.bw
	
	# Genome length
	genomelen=$(cat ./$chromtemp | awk '{total+=$2}END{print total}')
	
	# scale factor
	scale=$(awk 'BEGIN{printf "%.3f", '$genomelen'/'$totbase'}')
	echo -e "Here is the scale factor:" $scale
	
	# Write scaled bg
	awk '{printf "%s\t%s\t%s\t%.3f\n", $1, $2, $3, $4*'$scale'}' ./$plusbg > ./$norm_plusbg
	echo -e "Finished writing plus-strand bg file."
	awk '{printf "%s\t%s\t%s\t%.3f\n", $1, $2, $3, $4*'$scale'}' ./$minusbg > ./$norm_minusbg
	echo -e "Finished writing minus-strand bg file."
	
	# Convert to bw
	bedGraphToBigWig ./$norm_plusbg $chromtemp ./$norm_plusbw	
	echo -e "Finished writing plus-strand bw file."

	bedGraphToBigWig ./$norm_minusbg $chromtemp ./$norm_minusbw	
	echo -e "Finished writing minus-strand bw file."
	
	# Gzip
	gzip ./$norm_plusbg
	gzip ./$norm_minusbg
	
fi


# 9) Remove temp files
	rm ./$chromtemp ./$plustxpt ./$minustxpt ./$plusbg ./$minusbg ./$plustxptcov ./$minustxptcov


exit;

