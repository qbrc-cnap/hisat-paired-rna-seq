if(!require("DESeq2", character.only=T)) stop("Please install the DESeq2 package first.")

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
SAMPLE_ANNOTATION_FILE<-args[2]
CONDITION_A<-args[3]
CONDITION_B<-args[4]
OUTPUT_DESEQ_FILE <- args[5] 
OUTPUT_NORMALIZED_COUNTS_FILE <- args[6]
EDITED_ANNOTATIONS <- args[7]

# read the raw count matrix, genes as row names:
count_data <- read.table(RAW_COUNT_MATRIX, sep='\t', header = T, row.names = 1, stringsAsFactors = F)

# read the annotations
annotations <- read.table(SAMPLE_ANNOTATION_FILE, sep='\t', header = F, col.names = c('sample','condition'), stringsAsFactors = F)

# in case the samples had 'strange' names that R had to convert, alter the sample names in the annotation dataframe:
annotations[, 'sample'] <- make.names(annotations[, 'sample'])

# Keep only the rows that concern our contrast of interest
# Could keep it all together, but this is more explicit
selected_groups<-c(CONDITION_A, CONDITION_B)
annotations<-annotations[annotations$condition %in% selected_groups,]

# subset to keep only the samples in the count table.  This is important if the annotation
# file has more samples than the count matrix:
count_mtx_cols = colnames(count_data)
annotations <- annotations[annotations$sample %in% count_mtx_cols,]

# subset to only keep samples corresponding to the current groups in the count_data dataframe
count_data <- count_data[,annotations$sample]

# DESeq2 expects that the rownames of the annotation data frame are the sample names.  Set the rownames and drop that col
rownames(annotations) <- annotations$sample
annotations <- annotations[-1]

# Need to set the condition as a factor since it's going to be used as a design matrix
annotations$condition <- as.factor(annotations$condition)

# in case the names were changed by R, write a new sample annotation file so downstream processes don't get thrown off by the change:
write.table(annotations, file=EDITED_ANNOTATIONS, sep='\t', col.names=F, quote=F)

dds <- DESeqDataSetFromMatrix(countData = count_data,
							  colData = annotations,
							  design = ~condition)

dds <- DESeq(dds)
res <- results(dds, cooksCutoff=F, contrast=c("condition", CONDITION_B, CONDITION_A))
original_colnames = colnames(res)
n = length(original_colnames)
baseMeanA = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == CONDITION_A]) 
baseMeanB = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == CONDITION_B]) 
res = cbind(rownames(res), res[,1],baseMeanA, baseMeanB, as.data.frame(res[,2:n])) 
colnames(res) = c('Gene', 'overall_mean', CONDITION_A, CONDITION_B, original_colnames[2:n])
resOrdered <- res[order(res$padj),]
write.table(resOrdered, OUTPUT_DESEQ_FILE, sep='\t', quote=F, row.names=F)

#normalized counts:
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
nc <- cbind(gene=rownames(nc), nc)
write.table(nc, OUTPUT_NORMALIZED_COUNTS_FILE, sep='\t', quote=F, row.names=F)
