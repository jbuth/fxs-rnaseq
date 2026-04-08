## aggregate count matrix from salmon output from fxs-ms-rnaseq fastq files
# module load R/3.6.1
# R

# clear workspace
rm(list=ls()) 

# set the default for strings to be characters not factors
options(stringsAsFactors = FALSE)

# set directory as salmon output folder (each sample should have a folder here with quant.sf inside)
setwd("salmon")

## --- create tx2gene with biomart (transcript/gene annotations) --- ##

library(biomaRt)

bm = useEnsembl("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
a = listAttributes(bm); f= listFilters(bm);
tx2gene = getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"), mart=bm,useCache = FALSE)
dim(tx2gene) # [1] 144778      2 

files=paste0(list.files(),"/quant.sf") # sample_name/quant.sf
names(files)=list.files() # name each file so tximport adds names to count matrices

## --- import salmon output files to one matrix --- ##

# BiocManager::install("tximport")
library(tximport)

counts = tximport(
  files,                           # vector of filenames for transcript-level abundance
  type = c("salmon"),		           # options "salmon", "sailfish", "alevin", "kallisto", "rsem" etc
  txIn = TRUE,			               # TRUE means incoming file is transcript-level
  txOut = FALSE,			             # FALSE means will output to gene level, TRUE is both
  countsFromAbundance = c("no"),	 # options "no" (default), "scaledTPM", "lengthScaledTPM", etc
  tx2gene = tx2gene, 			         # two column data frame matching transcript.ids to gene.ids
  ignoreTxVersion = TRUE) 		     # ignoring version info

# counts is a list
# summary(counts)
    # Length  Class  Mode     
    # abundance           5867856 -none- numeric  
    # counts              5867856 -none- numeric  
    # length              5867856 -none- numeric  
    # countsFromAbundance       1 -none- character

# count matrix
# dim(counts$counts)
    # [1] 54332   108
    # 54,332 genes x 108 samples

## --- save count matrix and annotation --- ##

save(counts, tx2gene, file="data/fxs-rnaseq_aggregated_counts.RData")
