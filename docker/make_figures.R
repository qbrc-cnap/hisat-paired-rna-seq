#if(!require("reshape2", character.only=T)){install.packages("reshape2"); library(reshape2)}
#if(!require("ggdendro", character.only=T)){install.packages("ggdendro"); library(ggdendro)}
if(!require("gplots", character.only=T)){install.packages("gplots"); library(gplots)}

source("draw_heatmap.R")

# args from command line:
args<-commandArgs(TRUE)
DESEQ_OUTPUT_FILE<-args[1]
SAMPLE_ANNOTATION_FILE<-args[2]
NORM_COUNTS_FILE <- args[3]
OUTPUT_FIGURES_DIR <- args[4]
PADJ_THRESHOLD <- args[5]
LFC_THRESHOLD <- args[6]
CONTRAST_NAME <- args[7]
TOP_GENES_HM_SUFFIX <- args[8]
SIG_GENES_HM_SUFFIX <- args[9]

# Set threshold cutoff
thresholds=list(padj=as.numeric(PADJ_THRESHOLD),
  log2FoldChange=as.numeric(LFC_THRESHOLD))

# Set color Palette order
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")  # color blind friendly palette


# Read in results and store in results$dge, results$norm.mtx, and results$annotations
matrix_dt=read.table(NORM_COUNTS_FILE, header=T, row.names = 1, sep="\t")
dge_results=read.table(DESEQ_OUTPUT_FILE, header=T, sep="\t", stringsAsFactors = F)
annotations <- read.table(SAMPLE_ANNOTATION_FILE,
  sep='\t',
  header = F,
  col.names = c('Sample_ID','Group'), stringsAsFactors = F)

# convert any sample names that might have been changed by R:
annotations[, 'Sample_ID'] <- make.names(annotations[, 'Sample_ID'])

# have to remove any extra samples that might be in the annotation file:
annotations = annotations[annotations$Sample_ID %in% colnames(matrix_dt),]
rownames(annotations)=annotations$Sample

results=list(dge=dge_results,
             norm.mtx=matrix_dt,
             annotations=annotations)

# Make heatmap of top 20 genes
Top20.nData=results$norm.mtx[results$dge[1:40, "Gene"],]

## Add 1 to all values before log transformation
Top20.nData=log(Top20.nData+1)
Top20.nData<-as.matrix(Top20.nData)
## Draw Heatmap
top_heatmap_filepath = file.path(
  OUTPUT_FIGURES_DIR, 
  paste(CONTRAST_NAME, TOP_GENES_HM_SUFFIX, sep='.' )
)
Draw.Heatmap(Top20.nData,col.panel=c("blue", "white", "red"),
             scale="none",
             outfile=top_heatmap_filepath)

# Make heatmap for all significant genes
## Perform filtering based off pvalue
Sig.Gene=c(as.character(subset(results$dge,
                               padj<thresholds$padj &
                                 log2FoldChange>thresholds$log2FoldChange,
                               Gene)$Gene),
           as.character(subset(results$dge,
                               padj<thresholds$padj  &
                                 log2FoldChange< -thresholds$log2FoldChange,
                               Gene)$Gene))
sig.nData=as.data.frame(results$norm.mtx[Sig.Gene,])

sig_heatmap_filepath = file.path(
  OUTPUT_FIGURES_DIR, 
  paste(CONTRAST_NAME, SIG_GENES_HM_SUFFIX, sep='.' )
)

if(dim(sig.nData)[1] > 0){
  ## Add 1 to all values before log transformation
  sig.nData=log(sig.nData+1)
  sig.nData<-as.matrix(sig.nData)

  Draw.Heatmap(sig.nData, col.panel=c("blue", "white", "red"),
              scale="none",
              outfile=sig_heatmap_filepath)
} else {
    png(sig_heatmap_filepath)
    par(mar=c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x=0.5,y=0.5, sprintf('Zero genes passed significance\n at threshold of padj<%.2f, \nabs(log-fold-change) > %.2f',thresholds$padj,thresholds$log2FoldChange), cex = 1.6, col = "black")
    dev.off()
}


