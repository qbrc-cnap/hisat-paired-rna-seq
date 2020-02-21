workflow align_test {
    File r1_fastq
    File? r2_fastq
    File hisat_index_path
    File gtf
    String sample_name

    call perform_align {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            hisat_index_path = hisat_index_path,
            sample_name = sample_name
    }

    output {
        File bam = perform_align.sorted_bam
    }
}

task perform_align{
    # align to the reference genome using the STAR aligner
    # The STAR alignment produces a position-sorted BAM file

    # Input params passed by a parent Workflow:
    # We require that there is a single R1 fastq and possibly a single
    # R2 fastq for paired sequencing.
    # Genome is a string that matches one of the keys in the Map below
    # sample_name helps with naming files for eventual concatenation.
    File r1_fastq
    File? r2_fastq
    File hisat_index_path
    String sample_name

    # Default disk size in GB
    Int disk_size = 500

    command {
        set -euxo pipefail
        mkdir -p workspace/index
        #tar -xf ${hisat_index_path} -C workspace/index
        hisat2 [options]* -x workspace/index -1 ${r1_fastq} -2 ${r2_fastq} -S "${sample_name}.sam"
        samtools view -bh "${sample_name}.sam" > "${sample_name}.bam"
        samtools sort "${sample_name}.bam" > "${sample_name}.sort.bam"
    }

    output {
        File sorted_bam = "${sample_name}.sort.bam"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 8
        memory: "40 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

}
