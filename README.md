# Project: "pS22-LMNA"

This repository includes codes that Ikegami et al. used in the analysis presented in the paper.

## ChIP-seq

### Step 1: ChIPseq_pipeline_v2.sh 
Generates .pbs file to run in qsub system. This will take fastq, run bowtie, and produce bam.  

### Step 2: runMacs2.sh 
Takes bam and runs Tao's Macs2. Produces peakcalls (bed's) and fold-enrichment files (bg and bw). 

## RNA-seq

### Step 1: runTophat2_v2.sh
Takes fastq and produces bam.

### Step 2: makePerGeneSymbolTranscriptBed.sh
Takes ucsc table browser output and produces a file for one transcript per gene symbol. This file is often useful for chromatin-centric analyses of RNA.

### Step 3: computeRPKM_directional_v2.sh
Takes bam and transcript info and produces rpkm and base count.

## GRO-seq

### Step 1: bowtie_for_GROPE
Takes fastq and produces bam.

### Step 2: Bruseq_pipeline_v1.1_PE.sh
Takes bam and transcript info and produces normalized coverage per transcript and signal tracks for each strand.

## ATAC-seq

### Step 1: bowtie_for_ATAC
Takes fastq and produces bam. Uses first 38 nt.

### Step 2: atacbamToCutdens.sh
Takes sorted bam and returns Tn5 cutdensity profile for each replicate.

### Step 3: ATAC_peakcall
This decribes how to generate fold-enrichment bedgraph and call peaks. 




Kohta Ikegami
