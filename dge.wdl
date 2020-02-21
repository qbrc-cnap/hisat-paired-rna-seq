task run_differential_expression {

    File sample_annotations
    File raw_count_matrix
    String base_group
    String experimental_group
    String output_deseq2_suffix
    String normalized_counts_suffix
    String versus_sep
    Float padj_threshold
    Float lfc_threshold
    String top_genes_heatmap_suffix
    String sig_genes_heatmap_suffix

    Int disk_size = 30

    String contrast_name = experimental_group + versus_sep + base_group
    String output_deseq2 = contrast_name + "." + output_deseq2_suffix
    String normalized_counts = contrast_name + "." + normalized_counts_suffix
    String output_figures_dir = contrast_name + "_figures"
    String edited_annotations = 'edited_annotations.tsv'

    command <<<
        export OWD=$(pwd)
        mkdir ${output_figures_dir}

        Rscript /opt/software/deseq2.R \
            ${raw_count_matrix} \
            ${sample_annotations} \
            ${base_group} \
            ${experimental_group} \
            ${output_deseq2} \
            ${normalized_counts} \
            $OWD/${edited_annotations}

        # move the working dir to avoid headaches with R's source()
        cd /opt/software
        Rscript make_figures.R \
            $OWD/${output_deseq2} \
            $OWD/${edited_annotations} \
            $OWD/${normalized_counts} \
            $OWD/${output_figures_dir} \
            ${padj_threshold} \
            ${lfc_threshold} \
            ${contrast_name} \
            ${top_genes_heatmap_suffix} \
            ${sig_genes_heatmap_suffix}

        cd $OWD
        python3 /opt/software/make_dge_plots.py \
            -i ${output_deseq2} \
            -c ${normalized_counts} \
            -s $OWD/${edited_annotations} \
            -x ${contrast_name} \
            -o ${output_figures_dir} \
            -p ${padj_threshold} \
            -a ${base_group} \
            -b ${experimental_group}
    >>>

    output {
        File dge_table = "${output_deseq2}"
        File nc_table = "${normalized_counts}"
        Array[File] figures = glob("${output_figures_dir}/*")
    }  

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
