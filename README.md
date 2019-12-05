# pS22-LMNA

This repository includes codes that Ikegami et al. used in the analysis presented in the paper.

## ChIP-seq

### Step 1: ChIPseq_pipeline_v2.sh 
Generates .pbs file to run in qsub system. This will take fastq, run bowtie, and produce bam.  

### Step 2: runMacs2.sh 
Takes bam and runs Tao's Macs2. Produces peakcalls (bed's) and fold-enrichment files (bg and bw). 

## RNA-seq

### Step 1: runTophat2_v2.sh
Takes fastq and produces bam.




Kohta Ikegami
