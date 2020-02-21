task create_qc {
    #Array[File] star_logs
    Array[File] fc_logs
    Array[File] r1_fastqc_zips
    Array[File] dedup_metrics
    Array[File]? r2_fastqc_zips

    Int disk_size = 30

    command {
        multiqc .
    }

    output {
        File report = "multiqc_report.html"
    }
        
    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
