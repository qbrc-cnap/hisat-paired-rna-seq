Report for alignment and differential expression analysis
---

This document discusses the steps that were performed in the analysis pipeline.  It also describes the format of the output files and some brief interpretation.  For more detailed questions about interpretation of results, consult the documentation of the various tools.


## Results:

We summarize some brief results in this section.  Full results can be found in the files, as described in the *Outputs* section.  

**PCA** 

Principle component analysis (PCA) of the full normalized count matrix was performed. The first two components, PC1 and PC2, are shown.  Each sample groups are represented with a different color

![](other_figures/{{pca_plot}})

**Hierarchical clustering of samples**

Hierarchical clustering tree of all samples computed with euclidean distance considering all genes.

![](other_figures/{{hcl_plot}})

The following contrasts were performed, yielding the differentially expressed gene counts shown below.  The threshold for significance was set such that the adjusted p-value is less than {{adj_pval}}.  For the heatmap figures, the plotted genes were further limited to those with log-fold change magnitudes of {{lfc_threshold}} or greater.  When referenced, the "top" genes refers to the 40 genes with lowest p-value.   

|Experimental condition | Base condition| Upregulated | Downregulated | Result table | Heatmap of significant genes | Heatmap of top DE genes | Volcano |
|---|---|---|---|---|---|---|---|
{% for item in contrast_display %}
|{{item.exp_condition}} | {{item.base_condition}} | {{item.up_counts}}|{{item.down_counts}}| [Table](differential_expression/{{item.contrast_name}}.{{deseq2_output_file_suffix}}) | [Figure](differential_expression/{{item.contrast_name}}.{{sig_heatmap_file_suffix}}) | [Figure](differential_expression/{{item.contrast_name}}.{{top_heatmap_file_suffix}})|  [Figure](differential_expression/{{item.contrast_name}}.{{dynamic_volcano_file_suffix}})|
{% endfor %}

## Outputs:

This section describes the contents of the delivered results.

#### Alignments

Individual alignment files (in compressed BAM format, ending with "bam") are available for download, but are provided separately due to their typically large size.  If you download the BAM files, ensure that you also have the corresponding "index" files, which end with "bai".  Index files allow programs like IGV to use the BAM file in an efficient manner.

#### Main results

The main results are contained in a zip-archive and should be downloaded an "unzipped" on your local computer.  It contains several sub-directories which contain files produced in each step of the pipeline.

- **QC**
    - This directory contains an interactive HTML-based QC report which summarizes read quality, alignment quality, and other metrics.  It was produced by MultiQC, and information can be found at <https://multiqc.info/>.
    - Other QC plots are provided, produced by the RSeQC tool.  See documentation at <http://rseqc.sourceforge.net/> for details on each plot.
- **Quantifications**
    - Quantification tables, which give the number of reads aligned to each gene.  Files are tab-delimited.  These may be opened with your software of choice, including spreadsheet software such as Excel (note: <https://doi.org/10.1186/s13059-016-1044-7>).  For particulars on how this achieved, please see the featureCounts documentation or publication
- **Logs**
    - This contains logs and summaries produced by the various tools.  These can be used for troubleshooting, if necessary.
- **Differential expression results**

    Differential expression results are organized by the contrasts requested.  Thus, each folder, named by the contrast, contains the following files:

    - The differential expression results are summarized in a table in the file named `{{deseq2_output_file_suffix}}`.  It is saved in a tab-delimited text format.  Each row contains information about a particular gene and the results of the statistical test performed.  Note that we have added a couple of columns to the standard DESeq2 output.

    Although the particulars of DESeq2 are different, it can helpful to recall the Student's t-test when thinking about interpretation of the columns.  The parameter estimates and tests are different, but the ideas are very similar.  We are testing the null hypothesis that there is *no change* in the expression of a particular gene between the two conditions.

    The columns and brief interpretations are:

      - **Gene**: The gene symbol
      - **overall_mean**: You may think of this as a "blended mean" of the expression across all the samples considered (for this gene)
      - **Group1**: The average of the normalized expressions for samples from condition 1
      - **Group2**: The average of the normalized expressions for samples from condition 2
      - **log2FoldChange**: The logarithm (base 2) of the fold-change between the two sample groups.  Note that this is *not* simply based on the ratio of the average expressions in each group.  Details are in the DESEq2 publication.
      - **lfcSE**:  You may consider this as you would the standard error, which ultimately feeds into the p-value.  Consistent (i.e. not highly variable) expression within each group of samples yields lower standard error and makes it easier to determine if there are true differences in the mean expression of each group
      - **stat**: The value of the test statistic.
      - **pvalue**:  The raw p-value of the statistical test.  Lower values indicate more evidence to reject the null hypothesis.
      - **padj**: The "adjusted" p-value, which adjusts for the large number of statistical tests performed.  This addresses issues encountered in the "multiple-testing problem".  Corrections are based on the Benjamini-Hochberg procedure.

  - Normalized expression table

    A table of normalized read-counts (suitable for using when plotting expression) is provided in `{{normalized_counts_file_suffix}}`.  It is saved in a tab-delimited text file.  If multiple contrasts were performed, there will be differences in the tables of normalized counts, if different samples were used for each contrast.  However, we provide the normalized counts for each contrast since those were the expressions used in performing that particular differential contrast.     

- Figures

    We provide scatter plots of the top differentially expressed genes (if any) for quick reference.  Additionally, we provide a dynamic volcano plot which shows the log2 fold-change and adjusted p-value on a single plot.  This may be opened in any modern web browser.  Note that these figures are used for quick inspection, and it is not expected that such figures will be "publication-ready".  


## Methods:

Input fastq-format files are aligned to the {{genome}} reference genome using the HISAT2 aligner ({{hisat_version}}) [1].  BAM-format alignment files were filtered to retain only the primary-aligned reads using samtools ({{samtools_version}}) [3].  Additionally, "de-duplicated" versions of the primary-filtered BAM files were created using PicardTools' MarkDuplicates software ({{picard_mark_duplicates_version}})[4].  Both BAM files were indexed and quantified using featureCounts software ({{featurecounts_version}})[2] where counts were generated with respect to exon features.  Integer counts were concatenated into a file count "matrix" with rows denoting genes and samples denoting the samples.

Quality-control software included FastQC ({{fastqc_version}}), RSeQC ({{rseqc_version}}), and MultiQC ({{multiqc_version}}).  Please see the respective references for interpretation of output information and figures.

Note that we provide both the "unfiltered" and the "deduplicated" BAM files.  Depending on the quality of the experiment (e.g. very low input requiring many PCR cycles), it might makes sense to use the "deduplicated" version to reduce potential biases introduced by high-duplication rates.  By default, we only perform differential expression on expression counts derived from the "unfiltered" BAM files.

Integer read-count tables derived from RNA-seq alignments were analyzed for differential expression using Bioconductor's DESeq2 software.  Briefly, this software performs normalization to control for sequencing depth and subsequently performs differential expression testing based on a negative-binomial model.  For further details on both of these steps, please consult the documentation and publications for DESeq2 [5] and its older iteration, DESeq [6].

The R `sessionInfo()` produced the following output.  We print here so that the same combination of packages/software may be recreated, if necessary.

```
{{session_info}}
```

## Inputs:
The inputs to the workflow were given as:

The inputs to the workflow were given as:

Samples and sequencing fastq-format files:

{% for obj in file_display %}
  - {{obj.sample_name}}
    - R1 fastq: {{obj.r1}}
    - R2 fastq: {{obj.r2}}
{% endfor %}

Sample annotations file: `{{annotations_file}}`

Parsed sample and condition table:

|Sample|Condition|
|---|---|
{% for item in annotation_objs %}
|{{item.name}} | {{item.condition}} |
{% endfor %}


## Version control:
To facilitate reproducible analyses, the analysis pipeline used to process the data is kept under git-based version control.  The repository for this workflow is at 

<{{git_repo}}>

and the commit version was {{git_commit}}.

This allows us to run the *exact* same pipeline at any later time, discarding any updates or changes in the process that may have been added. 


#### References:

[1] Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019).
[2] Liao Y. and  Smyth G.K. and Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.  Bioinformatics. 2014

[3] Li H. and Handsaker B. and Wysoker A. and Fennell T. and Ruan J. and Homer N. and Marth G. and Abecasis G. and Durbin R. and 1000 Genome Project Data Processing Subgroup.  The Sequence alignment/map (SAM) format and SAMtools.  Bioinformatics. 2009.

[4] <http://broadinstitute.github.io/picard/>

[5] Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

[6] Anders S, Huber W (2010). "Differential expression analysis for sequence count data." Genome Biology, 11, R106. doi: 10.1186/gb-2010-11-10-r106.
