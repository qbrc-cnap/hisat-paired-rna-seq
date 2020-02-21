import "picard_tools.wdl" as picard_tools
import "hisat_align.wdl" as hisat_align
import "samtools.wdl" as samtools
import "feature_counts.wdl" as feature_counts
import "rseqc.wdl" as rseqc

workflow SingleSampleRnaSeqWorkflow{

    # This workflow operates on a single "sample" and
    # assumes that all the sequence reads are contained in 
    # 1 or 2 fastq files, for single- and paired-end sequencing,
    # respectively.
    #
    # Outputs are two count files, created from primary-filtered
    # and primary-filtered + deduplicated BAM files.

    File r1_fastq
    File? r2_fastq
    File hisat_index_path
    File gtf
    File bed_annotations

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1_fastq, "_R1.fastq.gz")

    # Perform the alignment, which outputs a sorted BAM
    call hisat_align.perform_align as alignment{
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            hisat_index_path = hisat_index_path,
            sample_name = sample_name
    }

    # Index that sorted BAM
    call samtools.samtools_index as index1 {
        input:
            input_bam = alignment.sorted_bam
    }

    # run the modified version of RseQC's infer experiment:
    #call rseqc.infer_experiment as infer_experiment{
    #    input:
    #        input_bam = alignment.sorted_bam,
    #        input_bam_index = index1.bam_index,
    #        bed_annotations = bed_annotations
    #}

    # run the remainder of the QC process
    call rseqc.qc_process as rseqc_process{
        input:
            input_bam = alignment.sorted_bam,
            input_bam_index = index1.bam_index,
    }

    # Filter for primary reads only
    call samtools.samtools_primary_filter as primary_filter{
        input:
            input_bam = alignment.sorted_bam,
            input_bam_index = index1.bam_index,
            sample_name = sample_name
    }

    # Index the primary-filtered BAM
    call samtools.samtools_index as index2 {
        input:
            input_bam = primary_filter.output_bam
    }

    # Mark and remove duplicates from the primary-filtered BAM
    call picard_tools.picard_deduplicate as deduplicate {
        input:
            input_bam = primary_filter.output_bam,
            input_bam_index = index2.bam_index
    }

    # Index the de-duplicated BAM
    call samtools.samtools_index as index3 {
        input:
            input_bam = deduplicate.output_bam
    }

    # Quantify the primary-filtered BAM
    call feature_counts.count_reads as quantify_primary {
        input:
            input_bam = primary_filter.output_bam,
            gtf = gtf,
            sample_name = sample_name,
            tag = "primary"
    }
    
    # Quantify the primary + deduplicated BAM
    call feature_counts.count_reads as quantify_deduplicated {
        input:
            input_bam = deduplicate.output_bam,
            gtf = gtf,
            sample_name = sample_name,
            tag = "primary_and_dedup"
    }

    output {
        File unfiltered_bam = alignment.sorted_bam
        File unfiltered_bam_index = index1.bam_index
        File primary_bam = primary_filter.output_bam
        File primary_bam_index = index2.bam_index 
        File primary_and_dedup_bam = deduplicate.output_bam
        File primary_and_dedup_bam_index = index3.bam_index
        File primary_filter_feature_counts_file = quantify_primary.count_output
        File primary_filter_feature_counts_summary = quantify_primary.count_output_summary
        File dedup_feature_counts_file = quantify_deduplicated.count_output
        File dedup_feature_counts_summary = quantify_deduplicated.count_output_summary
        File dedup_metrics = deduplicate.dedup_metrics
    }

}
