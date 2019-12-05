#!/bin/bash

# Version changes:
# 8-1-2016: Added a function to return the number of aligned reads (this version was temporally called v2.1)
# 12-05-2019: Removed email address

if [ $# -lt 10 ]; then # '-lt' stands for less than
	echo -e "\tby Kohta Ikegami 2016; all rights reserved\n"
	
	echo -e "\tWhat does this do?"
	echo -e "\tRuns bowtie2 on fastq files and spits MAPQ-filtered alignment (.bam), BedGraph (.bg), and Bigwig (.bw)\n"
	echo -e "\tUsage: `basename $0` [job name] [seed] [fastq] [genome] [ext (bp)] [MAPQ] [cpu] [memory (gb)] [time (hh)] [outdir (fullpath)]"
	echo -e "\tAll of the following 10 fields are required in order.\n"

	echo -e "\t[job name]\t The pbs file will be jobname.pbs"
	echo -e "\t[seed]\t\t Filename seed used in aligned bam, bedGraph, and bigwig. (e.g. KI310_H110_GM07492_hg19_BT2q40)"
	echo -e "\t[fastq]\t\t Accepts .gz. If multiple fastq, separate by commas."	
	echo -e "\t[genome]\t Genome to which reads are aligned (will look into /home/kikegami/ikgmk/computing_resources/genome_bfa/)"
	echo -e "\t[ext]\t\t Fragment size (bp) used to extend reads in BedGraph and Bigwig."
	echo -e "\t[MAPQ]\t\t Mapping quality. Use 40."
	echo -e "\t[cpu]\t\t Num of CPU cores. node => 4 processors (per node) => 16 cores (per processor)"
	echo -e "\t[memory]\t Memory usage in gb. e.g. 20 (=20gb)"
	echo -e "\t[time]\t\t Wall time by hour. e.g. 01 (=01:00:00)"
	echo -e "\t[outdir]\t Output directory in full path. No slash at the end.\n"

	exit 
fi

##################################
# 1. Set Variables
##################################

# 1) Variables

	genomebfa="/home/kikegami/ikgmk/computing_resources/genome_bfa"

	jobname=$1
	seed=$2			##e.g. KI310_H110_GM07492_hg19_BT2q40
	fastq=$3
	genome=$4
	ext=$5
	MAPQ=$6
	cpu=$7
	mem=$8"gb"
	runtime=${9}":00:00"
	outdir=${10}

	pbsfile=$jobname.pbs
	chromfile=$genomebfa/$genome"_chromLength_k1sort.txt"
	aligned=$seed.bam


#############################################
# 2. Check existig files and directories
#############################################

# 1) Outdir
	if [ ! -d $outdir ]; then mkdir $outdir; fi

#####################################################
# 3. Print Resource Manager Directives to .pbs file # 
#####################################################

echo '################# Resource Manager Directives #################'  > $outdir/$pbsfile  
echo ''

# Job name
echo '#PBS -N '$jobname'' >> $outdir/$pbsfile  

# Shell
echo '#PBS -S /bin/bash' >> $outdir/$pbsfile

# Email
# echo '#PBS -M email@uchicago.edu' >> $outdir/$pbsfile
# echo '#PBS -m ae' >> $outdir/$pbsfile 

# Run time
echo '#PBS -l walltime='$runtime'' >> $outdir/$pbsfile 

# Number of CPU core
echo '#PBS -l nodes=1:ppn='$cpu'' >> $outdir/$pbsfile 

# Memory usage
echo '#PBS -l mem='$mem'' >> $outdir/$pbsfile 

# Output
echo '#PBS -o '$outdir'' >> $outdir/$pbsfile
echo '#PBS -e '$outdir'' >> $outdir/$pbsfile

#####################################
# 4. Print commands to .pbs file # 
#####################################

echo '################# COMMANDS START HERE #################' >> $outdir/$pbsfile

# 0) Load modules
# Below is current setting using gardner
	echo 'module load gcc/6.2.0 bowtie2/2.1.0 bedtools/2.25.0 samtools/0.1.19 UCSCtools/343' >> $outdir/$pbsfile

# Below was setting for tarbell
#	echo 'module load bowtie2/2.1.0' >> $outdir/$pbsfile
#	echo 'module load bedtools/2.25.0' >> $outdir/$pbsfile
#	echo 'module load samtools/0.1.19' >> $outdir/$pbsfile
#	echo 'module load ucsc' >> $outdir/$pbsfile

# 1) Run Bowtie2
	
	## 1.1) Message
	echo 'echo -e "bowtie2 started at:"' >> $outdir/$pbsfile
	echo 'date' >> $outdir/$pbsfile
	echo 'echo -e "Here is the working directory:"' >> $outdir/$pbsfile
	echo 'pwd' >> $outdir/$pbsfile
	
	## 1.2) Bowtie2 run
	echo 'bowtie2 -p '$cpu' -x '$genomebfa'/'$genome' -U '$fastq' | samtools view -q '$MAPQ' -bS -o '$outdir'/'$aligned' -' >> $outdir/$pbsfile
	
	# Default alignment parameter is "--sensitive" (= -D 15 -R 2 -L 22 -i S,1,1.15)

	## 1.3) Get number of reads MAPQ-filtered aligned reads
	
	echo 'nreads=$(samtools view '$outdir'/'$aligned' -c)' >> $outdir/$pbsfile 
	echo 'echo -e "Here is the number of reads in '$aligned':" $nreads' >> $outdir/$pbsfile
	echo 'echo -e "ChIPseq_pipeline_v2.sh: Finished running Bowtie2."' >> $outdir/$pbsfile


# 2) Run bamToBgBw.sh

	## 2.1) Message
	echo 'echo -e "bamToBgBw started at:"' >> $outdir/$pbsfile
	echo 'date' >> $outdir/$pbsfile
	echo 'echo -e "Here is the working directory:"' >> $outdir/$pbsfile
	echo 'pwd' >> $outdir/$pbsfile	
	
	## 2.2) bamToBgBw run
	#	Usage: bamToBgBw.sh [bam] [file seed] [extension] [chromLength] [MAPQ] [outdir] [-c/n]
	echo 'bamToBgBw.sh '$outdir'/'$aligned' '$seed' '$ext' '$chromfile' '$MAPQ' '$outdir' -n' >> $outdir/$pbsfile

	echo 'echo -e "ChIPseq_pipeline_v2.sh: Finished running bamToBgBw.sh."' >> $outdir/$pbsfile

exit;

