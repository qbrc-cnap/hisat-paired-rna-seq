workflow test_rseqc{

    File input_bam
    File input_bam_index
    File bed_annotations

    call infer_experiment {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            bed_annotations = bed_annotations
    }

    output {
        File infer_experiment_result = infer_experiment.infer_results
    }

}

task infer_experiment {

    File input_bam
    File input_bam_index
    File bed_annotations

    Int disk_size = 100
    Int reads_sampled = 200000
    String outfile_name = "infer_experiment_output.csv"

    command {
        alternate_infer_experiment.py \
           -i ${input_bam} \
           -r ${bed_annotations} \
           -s ${reads_sampled} \
           -o ${outfile_name}
    }

    output {
        File infer_results = "${outfile_name}"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "8 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task qc_process {

    File input_bam
    File input_bam_index

    Int disk_size = 100

    command {
        echo "QC" > "qc_output.txt"
    }

    output {
        File qc_output = "qc_output.txt"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "8 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
