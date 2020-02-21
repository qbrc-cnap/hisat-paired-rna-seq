task generate_report {

    Array[String] r1_files
    Array[String] r2_files
    Array[File] deseq2_outputs
    File annotations
    String genome    
    String git_repo_url
    String git_commit_hash
    String normalized_counts_suffix
    String output_deseq2_suffix
    String versus_sep
    Float padj_threshold
    Float lfc_threshold
    String pca_filename
    String hc_tree_filename
    String top_genes_heatmap_suffix
    String sig_genes_heatmap_suffix

    String dynamic_volcano_file_suffix = "volcano.html"
    Int disk_size = 10

    command <<<

        # make a json file with various parameters:
        echo "{" >> config.json
        echo '"genome": "${genome}",' >>config.json
        echo '"pca_plot": "${pca_filename}",' >>config.json
        echo '"hcl_plot": "${hc_tree_filename}",' >>config.json
        echo '"adj_pval": "${padj_threshold}",' >>config.json
        echo '"lfc_threshold": "${lfc_threshold}",' >>config.json
        echo '"deseq2_output_file_suffix": "${output_deseq2_suffix}",' >>config.json
        echo '"normalized_counts_file_suffix": "${normalized_counts_suffix}",' >>config.json
        echo '"versus_sep": "${versus_sep}",' >>config.json
        echo '"git_repo": "${git_repo_url}",' >>config.json
        echo '"git_commit": "${git_commit_hash}",' >>config.json
        echo '"sig_heatmap_file_suffix": "${sig_genes_heatmap_suffix}",' >>config.json
        echo '"top_heatmap_file_suffix": "${top_genes_heatmap_suffix}",' >>config.json
        echo '"dynamic_volcano_file_suffix": "${dynamic_volcano_file_suffix}"}' >>config.json

        generate_report.py \
          -r1 ${sep=" " r1_files} \
          -r2 ${sep=" " r2_files} \
          -a ${annotations} \
          -d ${sep=" " deseq2_outputs} \
          -j config.json \
          -t /opt/report/report.md \
          -o completed_report.md

        pandoc -H /opt/report/report.css -s completed_report.md -o analysis_report.html
    >>>

    output {
        File report = "analysis_report.html"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
