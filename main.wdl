import "single_sample_rnaseq.wdl" as single_sample_rnaseq
import "feature_counts.wdl" as feature_counts
import "multiqc.wdl" as multiqc
import "fastqc.wdl" as fastqc
import "dge.wdl" as dge
import "figures.wdl" as figures
import "report.wdl" as reporting


workflow HisatPairedRnaSeqAndDgeWorkflow{
    # This workflow is a 'super' workflow that parallelizes
    # RNA-seq analysis over multiple samples

    Array[File] r1_files
    Array[File] r2_files
    File sample_annotations
    Array[String] base_conditions
    Array[String] experimental_conditions
    String genome
    File hisat_index_path
    File gtf
    File bed_annotations
    String output_zip_name
    String git_repo_url
    String git_commit_hash

    Float padj_threshold = 0.01
    Float lfc_threshold = 1.5

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    Array[Pair[String, String]] contrast_pairs = zip(base_conditions, experimental_conditions)

    String versus_sep = "_versus_"
    String normalized_counts_suffix = "normalized_counts.tsv"
    String output_deseq2_suffix = "deseq2_results.tsv"
    String pca_filename = 'pca.png'
    String hc_tree_filename = 'hctree.png'
    String top_genes_heatmap_suffix = 'top_genes_heatmap.png'
    String sig_genes_heatmap_suffix = 'significant_genes_heatmap.png'

    call reformat_annotations {
        input:
            sample_annotations = sample_annotations
    }

    scatter(item in fastq_pairs){

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
        }

        call single_sample_rnaseq.SingleSampleRnaSeqWorkflow as single_sample_process{
            input:
                r1_fastq = item.left,
                r2_fastq = item.right,
                hisat_index_path = hisat_index_path,
                gtf = gtf,
                bed_annotations = bed_annotations
        }
    }

    call feature_counts.concatenate as merge_primary_counts {
        input:
            count_files = single_sample_process.primary_filter_feature_counts_file,
            output_filename = "raw_primary_counts.tsv"
    }

    call feature_counts.concatenate as merge_dedup_counts {
        input:
            count_files = single_sample_process.dedup_feature_counts_file,
            output_filename = "raw_primary_and_deduplicated_counts.tsv"
    }

    call multiqc.create_qc as experimental_qc {
        input:
            #star_logs = single_sample_process.star_log,
            fc_logs = single_sample_process.primary_filter_feature_counts_summary,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            dedup_metrics = single_sample_process.dedup_metrics,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    call figures.create_contrast_independent_figures as make_figs {
        input:
            sample_annotations = reformat_annotations.reformatted_annotations,
            raw_count_matrix = merge_primary_counts.count_matrix,
            pca_filename = pca_filename,
            hc_tree_filename = hc_tree_filename      
    }

    scatter(item in contrast_pairs){
        call dge.run_differential_expression as run_dge {
            input:
                sample_annotations = reformat_annotations.reformatted_annotations,
                raw_count_matrix = merge_primary_counts.count_matrix,
                base_group = item.left,
                experimental_group = item.right,
                output_deseq2_suffix = output_deseq2_suffix,
                normalized_counts_suffix = normalized_counts_suffix,
                versus_sep = versus_sep,
                padj_threshold = padj_threshold,
                lfc_threshold = lfc_threshold,
                top_genes_heatmap_suffix = top_genes_heatmap_suffix,
                sig_genes_heatmap_suffix = sig_genes_heatmap_suffix
        }
    }

    call reporting.generate_report as generate_report{
        input:
            r1_files = r1_files,
            r2_files = r2_files,
            genome = genome,
            git_commit_hash = git_commit_hash,
            git_repo_url = git_repo_url,
            annotations = reformat_annotations.reformatted_annotations,
            deseq2_outputs = run_dge.dge_table,
            output_deseq2_suffix = output_deseq2_suffix,
            normalized_counts_suffix = normalized_counts_suffix,
            versus_sep = versus_sep,
            padj_threshold = padj_threshold,
            lfc_threshold = lfc_threshold,
            pca_filename = pca_filename,
            hc_tree_filename = hc_tree_filename,
            sig_genes_heatmap_suffix = sig_genes_heatmap_suffix,
            top_genes_heatmap_suffix = top_genes_heatmap_suffix
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            primary_fc_file = merge_primary_counts.count_matrix,
            dedup_fc_file = merge_dedup_counts.count_matrix,
            primary_bam_files = single_sample_process.primary_bam,
            primary_bam_index_files = single_sample_process.primary_bam_index,
            #star_logs = single_sample_process.star_log,
            dedup_fc_summaries = single_sample_process.dedup_feature_counts_summary, 
            primary_fc_summaries = single_sample_process.primary_filter_feature_counts_summary,
            dedup_metrics = single_sample_process.dedup_metrics,
            multiqc_report = experimental_qc.report,
            analysis_report = generate_report.report,
            deseq2_outputs = run_dge.dge_table,
            normalized_counts_files = run_dge.nc_table,
            contrast_figures = run_dge.figures,
            pca_figure = make_figs.pca,
            hctree_figure = make_figs.hctree
    }

    output {
        File zip_out = zip_results.zip_out
        File alignments_zip = zip_results.alignments_zip
    }

    meta {
        workflow_title : "Paired-end RNA-Seq basic differential expression"
        workflow_short_description : "For determining differential expression from a basic paired-end RNA-seq experiment"
        workflow_long_description : "Use this workflow for aligning with HISAT2, quantifying, and testing differential expression with DESeq2 from a paired-end RNA-seq experiment."
    }
}

task reformat_annotations {

    File sample_annotations
    String output_file = "reformatted_annotations.tsv"
    Int disk_size = 10
    
    command {
        python3 /opt/software/reformat_annotations.py -a ${sample_annotations} -o ${output_file}
    }

    output {
        File reformatted_annotations = "${output_file}"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "2 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}


task zip_results {

    String zip_name 

    File primary_fc_file
    File dedup_fc_file
    Array[File] primary_bam_files
    Array[File] primary_bam_index_files
    Array[File] hisat_logs
    Array[File] dedup_fc_summaries 
    Array[File] primary_fc_summaries
    Array[File] dedup_metrics
    File multiqc_report
    File analysis_report
    Array[File] deseq2_outputs
    Array[File] normalized_counts_files
    Array[Array[File]] contrast_figures
    File pca_figure
    File hctree_figure

    Array[File] contrast_figure_list = flatten(contrast_figures)

    Int disk_size = 1000

    command {

        mkdir alignments
        mv -t alignments ${sep=" " primary_bam_files}
        mv -t alignments ${sep=" " primary_bam_index_files}
        zip -r "alignments.zip" alignments

        mkdir report
        mkdir report/quantifications
        mkdir report/qc
        mkdir report/logs
        mkdir report/differential_expression
        mkdir report/other_figures

        mv ${primary_fc_file} report/quantifications/
        mv ${dedup_fc_file} report/quantifications/
        mv ${multiqc_report} report/qc/
        mv -t report/logs ${sep=" " hisat_logs}
        mv -t report/logs ${sep=" " dedup_fc_summaries}
        mv -t report/logs ${sep=" " primary_fc_summaries}
        mv -t report/logs ${sep=" " dedup_metrics}
        mv -t report/differential_expression ${sep=" " deseq2_outputs}
        mv -t report/differential_expression ${sep=" " normalized_counts_files}
        mv -t report/differential_expression ${sep=" " contrast_figure_list}

        mv ${pca_figure} report/other_figures/
        mv ${hctree_figure} report/other_figures/

        mv ${analysis_report} report/
        zip -r "${zip_name}.zip" report
    }

    output {
        File zip_out = "${zip_name}.zip"
        File alignments_zip = "alignments.zip"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
