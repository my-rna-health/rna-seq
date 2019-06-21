version development

workflow StarIndex {

    input {
        File referenceGenome
        String indexDir
        File gtf
        Int threads = 8
    }

    call star_index {
        input:
            path = indexDir,
            genomeFasta = referenceGenome,
            threads = threads,
            gtf = gtf
    }

    call copy_folder_content {
        input:
            folder = star_index.folder,
            destination = indexDir
    }


    output {
        Directory out =  copy_folder_content.out
    }
}

task star_index {
    input {
        String path
        File genomeFasta
        File gtf
        Int threads
        String ram = "58000000000"
        Int bins = 16
    }

    command {
        mkdir -p ~{path}
        /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --runMode genomeGenerate \
        --genomeDir ~{path} \
        --genomeFastaFiles ~{genomeFasta}  \
        --sjdbGTFfile ~{gtf} \
        --limitGenomeGenerateRAM ~{ram} \
        --genomeChrBinNbits ~{bins}
    }

    runtime {
        docker: "quay.io/biocontainers/star@sha256:6556f4d3a9f767f93cd8841ef0ac0f15d29d2c69b6a6f1f9b6ccc3e1da40207b" #2.7.1a--0
      }


    output {
        Directory folder = path
    }

}

task copy_folder_content {
    input {
        Directory folder
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{folder}/. ~{where}
    }

    output {
        Directory out = where
    }
}