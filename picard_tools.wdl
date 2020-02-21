workflow picard_test {
    File input_bam
    File input_bam_index

    call picard_deduplicate {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index
    }

    output {
        File output_bam = picard_deduplicate.output_bam
        File metrics = picard_deduplicate.dedup_metrics
    }
}

task picard_deduplicate {

    File input_bam
    File input_bam_index

    # construct the name by appending onto existing:
    String output_bam_basename = basename(input_bam)
    String output_bam_name = sub(output_bam_basename, "\\.bam", ".duplicates_removed.bam")

    Int disk_size = 300

    command {
        mkdir tmp
        export TMPDIR=$(pwd)/tmp
        java -Xmx28g -jar $PICARD_JAR MarkDuplicates \
	      INPUT=${input_bam} \
	      OUTPUT="${output_bam_name}" \
	      ASSUME_SORTED=TRUE \
	      TMP_DIR=$TMPDIR \
	      REMOVE_DUPLICATES=TRUE \
	      METRICS_FILE="${output_bam_name}".metrics.out \
          VALIDATION_STRINGENCY=LENIENT
    }

    output {
        File output_bam = "${output_bam_name}"
        File dedup_metrics = "${output_bam_name}.metrics.out"
    }

    runtime {
        docker: "docker.io/blawney/hisat_rnaseq:v0.0.1"
        cpu: 8
        memory: "32 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

}
