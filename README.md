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

## RNA-seq downstream processing

### Combat run

```r
# 1) Compute min 

my.rowmin <- aging.rpkm %>% 
	dplyr::select(s42:KI436) %>%
	t() %>% as.data.frame %>%
	summarise_all(funs(min)) %>%
	t()
 
# 2) Qnorm rpkm
 
my.subset <- aging.rpkm %>% 
	dplyr::bind_cols(data.frame(rpkm.all.min=my.rowmin)) %>%
	dplyr::left_join(rpkm4.genc1 %>% dplyr::select(transcript.id, gene.type)) %>%
	dplyr::filter(gene.type == "protein_coding", rpkm.all.min > 0.01) %>%
	dplyr::select(-s78, -KI429) %>%
	dplyr::select(
		transcript.id,
		s75:s84,
		s42:s51,
		KI430,KI433,KI434,
		KI431,KI432,KI435,KI436
    ) %>%
    dplyr::mutate_at(vars(s75:KI436), funs(log10(.)))

# > dim(my.subset)
# [1] 11613    27

my.subset.q <- my.subset %>% 
	dplyr::select(-transcript.id) %>%
	as.matrix %>% 
	normalize.quantiles %>%
	as.data.frame
	
colnames(my.subset.q) <- names(my.subset %>% dplyr::select(-transcript.id))
rownames(my.subset.q) <- my.subset$transcript.id
 
# 2) colData 
colData <- data.frame(
    cell.type=c( rep("Normal",9), rep("HGPS",10), rep("Normal",3), rep("HGPS",4) ),
    study.type=c( rep("aging",19), rep("ki",7) ),
    row.names= my.subset.q %>% names
)

my.batch = colData$study.type
my.modcombat = model.matrix(~1, data=colData)

# 4) Combat
my.subset.q.combat = ComBat(
	dat=my.subset.q, 
	batch=my.batch, 
	mod=NULL, 
	par.prior=TRUE, 
	prior.plots=TRUE
)

```

## DESeq2

```r
# > package.version("DESeq2")
# [1] "1.14.1"

# 1) Remove two outliers
my.subset <- aging.count %>%
	dplyr::select(-s78, -KI429) %>%
	dplyr::select(
		transcript.id,
		s75:s84,
		s42:s51,
		KI430,KI433,KI434,
		KI431,KI432,KI435,KI436
	)

#> names(my.subset)
# [1] "transcript.id" "s75"           "s76"           "s77"           "s79"           "s80"           "s81"          
# [8] "s82"           "s83"           "s84"           "s42"           "s43"           "s44"           "s45"          
# [15] "s46"           "s47"           "s48"           "s49"           "s50"           "s51"           "KI430"        
# [22] "KI433"         "KI434"         "KI431"         "KI432"         "KI435"         "KI436"        

# 2) colData 
colData <- data.frame(
	cell.type=c( rep("Normal",9), rep("HGPS",10), rep("Normal",3), rep("HGPS",4) ),
	study.type=c( rep("aging",19), rep("ki",7) ),
	row.names= my.subset %>% dplyr::select(-transcript.id) %>% names
)

# > colData
#       cell.type study.type
# s75      Normal      aging
# s76      Normal      aging
# s77      Normal      aging
# s79      Normal      aging
# s80      Normal      aging
# s81      Normal      aging
# s82      Normal      aging
# s83      Normal      aging
# s84      Normal      aging
# s42        HGPS      aging
# s43        HGPS      aging
# s44        HGPS      aging
# s45        HGPS      aging
# s46        HGPS      aging
# s47        HGPS      aging
# s48        HGPS      aging
# s49        HGPS      aging
# s50        HGPS      aging
# s51        HGPS      aging
# KI430    Normal         ki
# KI433    Normal         ki
# KI434    Normal         ki
# KI431      HGPS         ki
# KI432      HGPS         ki
# KI435      HGPS         ki
# KI436      HGPS         ki

# 3) Data formatting
rna.norm.hgps.de <- DESeqDataSetFromMatrix(
        countData=data.frame(
            my.subset %>% dplyr::select(-transcript.id),
            row.names=my.subset$transcript.id
        ),      
        colData=colData, 
        design= ~ study.type + cell.type
)

# 4) DESeq run
rna.norm.hgps.de <- DESeq(rna.norm.hgps.de)

# out of 22966 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 1521, 6.6% 
# LFC < 0 (down)   : 1418, 6.2% 
# outliers [1]     : 1961, 8.5% 
# low counts [2]   : 7341, 32% 
# (mean count < 326)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# 5) Norm vs HGPS
rna.norm.hgps.de05 <- extract.deseq(rna.norm.hgps.de, "rna.norm.hgps.de05", 0.05, "cell.type", "HGPS", "Normal")

########################################################
# extract.deseq
# Function to extract deseq data for later dplyr use
########################################################

extract.deseq <- function(dds,seed,alpha,condition,contrast1,contrast2) {
	
	# Example parameters and variables
	# dds: eg1.wt.ko.de
	# seed: "eg1.wt.ko.de"
	# alpha: 0.05
	# condition: "celltype"
	# contrast1: "KO"
	# contrast2: "WT" 

	library(DESeq2);
	
	x <- results(dds, alpha=alpha, contrast=c(condition, contrast1, contrast2));

	summary(x);
	
	x <- data.frame(x) %>%
	dplyr::mutate(transcript.id=factor(rownames(x))) %>%
	dplyr::select(
		transcript.id,
		lfc=log2FoldChange,
		padj=padj
	) %>%
	dplyr::mutate(bin=
		as.factor(
			ifelse(is.na(lfc) | is.na(padj), "no",
				ifelse(padj < alpha & lfc > 0, "up",
					ifelse(padj < alpha & lfc < 0, "dn", "no")
				)
			)
		)
	) %>%
	dplyr::mutate(bin=factor(bin, levels=c("up","dn","no")))
	
	names(x) <- c("transcript.id", paste0(seed,".lfc"), paste0(seed,".padj"), paste0(seed,".3bin"))

	return(x);
}
```




## Heatmap on differentially expressed genes


```r
# 1) Dataset

my.subset <- aging.rpkm %>%
	dplyr::filter(combat.2bin=="yes") %>%
	dplyr::filter(
		rna.norm.hgps.de05.lfc05.3bin != "no"
	) 	

# > dim(my.subset)
# [1] 1117   88

# 2) Plotdata

my.plotdat <- my.subset %>%
	dplyr::select(ends_with("_combat")) %>%
	dplyr::rename(
		AG11513.s42=s42_combat,
		HGADFN167.s43=s43_combat,
		HGADFN188.s44=s44_combat,
		HGADFN127.s45=s45_combat,
		HGADFN164.s46=s46_combat,
		HGADFN169.s47=s47_combat,
		HGADFN178.s48=s48_combat,
		HGADFN122.s49=s49_combat,
		HGADFN143.s50=s50_combat,
		HGADFN367.s51=s51_combat,
		GM00969.s75=s75_combat,
		GM05565.s76=s76_combat,
		GM00498.s77=s77_combat,
		GM05400.s79=s79_combat,
		GM05757.s80=s80_combat,
		GM00409.s81=s81_combat,
		GM00499.s82=s82_combat,
		GM08398.s83=s83_combat,
		GM00038.s84=s84_combat,
		GM07492.KI430=KI430_combat,
		AG11498.KI431=KI431_combat,
		HGADFN167.KI432=KI432_combat,
		GM08398.KI433=KI433_combat,
		GM07492.KI434=KI434_combat,
		AG11498.KI435=KI435_combat,
		HGADFN167.KI436=KI436_combat
	) %>%
	t %>% scale() %>% t    
row.names(my.plotdat) <- row.names(my.subset)

# 3) Annotations
annotation_col <- data.frame(
		cell.type=c( rep("Normal",9), rep("HGPS",10), rep("Normal",3), rep("HGPS",4) ),
		study.type=c( rep("aging",19), rep("ki",7) ),
		row.names= colnames(my.plotdat)
	)

annotation_row <- my.subset %>% 
		dplyr::select(rna.norm.hgps.de05.lfc05.3bin, rna.norm.hgps.de05.lfc05.gm07ag11.3bin)

# 4) Annotation colors
col4 <-  c(brewer.pal(n = 11, name = "Spectral")[c(3,10)],"gray80","black")
bicol <- hue_pal()(2) 
study.col <- brewer.pal(n = 8, name = "Dark2") 

annotation_color <- list(
		rna.norm.hgps.de05.lfc05.3bin = c(up=col4[1], dn=col4[2], no=col4[3]),
		study.type= c(ki=study.col[5], aging=study.col[6]),
		cell.type= c(Normal=bicol[2], HGPS=bicol[1])
	)

# Plot
pheatmap(my.plotdat,
	color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
	breaks=seq(-2, 2, length.out = 101),
	clustering_distance_cols = "correlation",
	annotation_col = annotation_col,
	annotation_row = annotation_row,
	annotation_colors = annotation_color,
	show_rownames =F,
	border_color = NA,
	cellwidth=10,
	cellheight=0.3,
	filename="heatmap.pdf",
	silent =T,
	main="rna.norm.hgps.de05.lfc05.3bin != no\ncombat2bin=yes"
)
dev.off()

```




Kohta Ikegami
