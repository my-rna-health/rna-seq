version development

workflow methyldackel {
    input {
        File genome
        File bam
        Int threads = 4
    }

    call methyldackel {
        input:
            bam = bam,
            genome = genome,
            threads = threads
    }
}

task methyldackel {

    input {
        File bam
        File genome
        Int threads = 4
    }

    command {
        MethylDackel extract --CHH --CHG --counts -@ ~{threads} ~{genome} ~{bam}
    }

    runtime {
        docker: "quay.io/biocontainers/methyldackel@sha256:d434c3e320a40648a3c74e268f410c57649ab208fcde4da93677243b22900c55" #0.3.0--h84994c4_3
    }

    output {
        File cpg = "alignments_CpG.bedGraph"
        File counts = "alignments.counts.bedGraph"
        #File chh = "-CHH and --CHG
    }
}