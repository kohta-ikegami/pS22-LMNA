#!/bin/bash

## Usage: makePerGeneSymbolTranscriptBed.sh [ucsc_download.txt] > [stdout]
## Requirement: bedtools v2.25.0

if [ $# -lt 1 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2016"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tGenerate a list of transcripts representing each gene symbol.\n"
	echo -e "Usage:  \t`basename $0` [ucsc_download.txt] [spike-in.gtf] > [stdout]"
	echo -e "        \t- [spike-in.gtf] is optional."
	echo -e "\n"
	exit 
fi

# Instruction on downloading GencodeBasicV19 file from UCSC

	# 1) Go to: http://genome.ucsc.edu/cgi-bin/hgTables
	# 2) Select fields as below:
		# clade: Mammal
		# genome: Human
		# assembly: Feb. 2009 (GRCh37/hg19)
		# group: Genes and Gene Predictions
		# track: GENCODE Genes V19
		# table: Basic (wgEncodeGencodeBasicV19)
		# region: genome
		# identifiers: (not selected)
		# filtered: (not selected)
		# subtrack merge: (not selected)
		# intersection: (not selected)
		# correlation: (not selected)
		# output format: selected fields from primary and related tables
		# output file: [filename to be saved]
		# file type returned: plain text
	# 3) click [get output]
	# 4) In the next page, 
		# under "Select Fields from hg19.wgEncodeGencodeBasicV19", check
			# name, chrom, strand, txStart, exEnd, exonCount, exonStarts, exonEnds, name2
		# "hg19.wgEncodeGencodeAttrsV19 fields", check: 
			# geneType, transcriptType
	# 5) Click [get output]

# Files
	tablebrowserfile=$1

# Processing

	# 0) Filed names in the downloaded file
		# - 1. hg19.wgEncodeGencodeBasicV19.name		(<== transcript.id)
		# - 2. hg19.wgEncodeGencodeBasicV19.chrom	
		# - 3. hg19.wgEncodeGencodeBasicV19.strand
		# - 4. hg19.wgEncodeGencodeBasicV19.txStart
		# - 5. hg19.wgEncodeGencodeBasicV19.txEnd
		# - 6. hg19.wgEncodeGencodeBasicV19.exonCount
		# - 7. hg19.wgEncodeGencodeBasicV19.exonStarts
		# - 8. hg19.wgEncodeGencodeBasicV19.exonEnds
		# - 9. hg19.wgEncodeGencodeBasicV19.name2		(<== gene.symbol)	
		# - 10. hg19.wgEncodeGencodeAttrsV19.geneType
		# - 11. hg19.wgEncodeGencodeAttrsV19.transcriptType

	# 1) Only include transcripts whose biotype matches that of the gene (99901 => 87336).
	awk '$10==$11' $tablebrowserfile |
	
	# 2) Only include transcripts whose (gene's) biotypes are antisense, or protein_coding, or lincRNA (87336 => 76094).
	awk '$10=="antisense" || $10=="protein_coding" || $10=="lincRNA"' |
	
	# 3) Only include transcripts whose transcript_id appear only once in the list (76094 => 75968). 
	sort -f -k1,1 | bedtools groupby -g 1 -c 1 -o count -full | awk '$12 == 1' | cut -f1-11 |
	
	# 4) Compute transcript length
	awk '{OFS="\t"; print $0,$5-$4}' |
	
	# 5) Sort by gene.symbol, exonCount (largest first), length (largest first), transcript.id in that order (smallest in the alphanumerical order first)
	sort -f -k9,9 -k6,6nr -k12,12nr -k1,1 |
	
	# 6) For each gene.symbol, choose the first line. This transcript has the largest exonCount, and if it ties, the largest length, and if it ties, the smallest transcript.id in the alphanumeric order (75968 => 32517). 
	bedtools groupby -g 9 -c 9 -o first -full |
	
	# 7) Output as a temp file with field names (32517 => 32518).
	awk '{OFS="\t"; print $2,$4,$5,$1,$12,$3,$9,$10,$6,$7,$8}' > PerGeneSymbolTranscript.temp

	# 8) Add ERCC92 spike-in if included.
	if [ -n "$2" ] ; then
		awk -F ";|\t|\"" '{OFS="\t"; print $1, $4, $5, $13, $5-$4, $7, $10, $2, 1, $4, $5}' $2 |
		cat - PerGeneSymbolTranscript.temp
	else
		cat PerGeneSymbolTranscript.temp 
	fi |
	
	# 9) Sort by transcript.id and output
	sort -f -k4,4 |
	awk 'BEGIN{print "chrom\ttxStart\ttxEnd\ttranscript.id\tlength\tstrand\tgene.symbol\tgeneType\texonCount\texonStarts\texonEnds"}{OFS="\t"; print $0}' 
	
	rm PerGeneSymbolTranscript.temp

exit

