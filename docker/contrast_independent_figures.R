#if(!require("reshape2", character.only=T)){install.packages("reshape2"); library(reshape2)}
if(!require("ggdendro", character.only=T)){install.packages("ggdendro"); library(ggdendro)}
if(!require("gplots", character.only=T)){install.packages("gplots"); library(gplots)}

library(ggplot2)
library(DESeq2)
source("draw_pca.R")

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
SAMPLE_ANNOTATION_FILE<-args[2]
PCA_FILENAME <- args[3]
HC_TREE_FILENAME <- args[4]

# Need the normalized counts:
count_data <- read.table(RAW_COUNT_MATRIX,
    sep='\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F)

# read the annotations
annotations <- read.table(SAMPLE_ANNOTATION_FILE,
    sep='\t',
    header = F,
    col.names = c('Sample_ID','Group'),
    stringsAsFactors = F)

# in case the names were not legal for R, change the annotations 
annotations[,'Sample_ID'] = make.names(annotations[,'Sample_ID'])

# subset to keep only the samples in the count table.  This is important if the annotation
# file has more samples than the count matrix:
count_mtx_cols = colnames(count_data)
annotations <- annotations[annotations$Sample_ID %in% count_mtx_cols,]

# need to harmonize the order of rows/cols:
count_data <- count_data[,annotations$Sample_ID]

# DESeq2 expects that the rownames of the annotation data frame are the sample names.
rownames(annotations) <- annotations$Sample_ID

# Need to set the condition as a factor since it's going to be used as a design matrix
annotations$Group <- as.factor(annotations$Group)

# since we are not actually performing differential expression, give
# the design as simply ~1
dds <- DESeqDataSetFromMatrix(countData = count_data, colData=annotations, design=~1)
dds = estimateSizeFactors(dds)
nc = counts(dds, normalized=T)

# Set color Palette order
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")  # color blind friendly palette

results=list( norm.mtx=nc,
             annotations=annotations)

# Make PCA plot and store in figures/pca.png
pca_plot=Draw.PCA(results=results, res=600, width=6, 
                  height=4, group.colors=cbPalette,
                  outfile=PCA_FILENAME)

# Make a hctree on the full matrix
hctree=as.dendrogram(hclust(dist(t(results$norm.mtx), method='euclidean')))
hctree_dend=ggdendrogram(hctree)
ggsave(HC_TREE_FILENAME, plot=hctree_dend, dpi=300 )
