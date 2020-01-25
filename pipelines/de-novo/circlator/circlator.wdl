version development

workflow circlator {
    input{
        Array[File]+ reads
        File assembly
        String out
    }

    call circle{
        input: reads = reads,
        assembly = assembly,
        directory = out
    }

    call copy{
        input
    }

    output {
        Directory out = out
    }
}



task circle {

    input {
        Array[File]+ reads
        File assembly
        String directory
    }


    command {
        circlator all ~{assembly} ~{reads} ~{directory}
    }

    runtime {
        docker: "quay.io/biocontainers/circlator@sha256:e1fb7edcd5323e82c765bd2cc30085e6af51d397d555037d8d3948f1dd33af84"
    }
}