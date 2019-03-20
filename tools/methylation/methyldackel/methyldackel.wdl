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

    String name = "alignment"

    command {
        MethylDackel extract --CHH --CHG --counts -@ ~{threads} ~{genome} ~{bam}
    }

    runtime {
        docker: "quay.io/biocontainers/methyldackel@sha256:d434c3e320a40648a3c74e268f410c57649ab208fcde4da93677243b22900c55" #0.3.0--h84994c4_3
    }

    output {
        File chg = name + "_CHG.counts.bedGraph"
        File chh = name + "_CHH.counts.bedGraph"
        File cpg = name + "_CpG.counts.bedGraph"
    }
}