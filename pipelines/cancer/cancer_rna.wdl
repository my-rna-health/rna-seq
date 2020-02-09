version development

workflow cancer_rna{
    input {
        Array[File] reads
        Directory index
    }

    call star_fusion{
        input:
            left = reads[0],
            right = reads[1],
            index = index

    }
}

task star_fusion {
    input {
        File left
        File right
        Directory index
    }

    command {
        /usr/local/src/STAR-Fusion/STAR-Fusion \
            --left_fq ~{left} \
            --right_fq ~{right} \
            --genome_lib_dir ~{index} \
            -O output \
            --FusionInspector validate \
            --examine_coding_effect \
            --denovo_reconstruct
    }

    runtime {
        docker: "trinityctat/ctatfusion:latest"
    }

    output {
        Directory out = "output"
    }
}