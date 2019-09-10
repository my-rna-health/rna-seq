version development

workflow StarIndex {

    input {
        File referenceGenome
        String indexDir
        File gtf
        Int threads = 4
        String max_ram = "50000000000"
        Int bins = 15
    }

    call star_index {
        input:
            path = indexDir,
            genomeFasta = referenceGenome,
            threads = threads,
            gtf = gtf,
            ram = max_ram,
            bins = bins
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
        String ram
        Int bins
    }

    String dirname = basename(path)

    command {
        mkdir ~{dirname}
        /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --runMode genomeGenerate \
        --genomeDir ~{dirname} \
        --genomeFastaFiles ~{genomeFasta}  \
        --sjdbGTFfile ~{gtf} \
        --limitGenomeGenerateRAM ~{ram} \
        --genomeChrBinNbits ~{bins}
    }

    runtime {
         docker: "quay.io/biocontainers/star@sha256:f9b0406354ff2e5ccfadaef6fde6367c7bcb4bdc7e67920f0f827a6ff6bf4fb5" #2.7.2b--0
         memory_mb: "51200"
      }


    output {
        Directory folder = dirname
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