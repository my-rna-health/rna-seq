workflow StarIndex {
    File referenceGenome
    File indexDir
    Int threads
    Int binBits

    call star_index {
        input:
            genomeDir = indexDir,
            genomeFasta = referenceGenome,
            threads = 8,
            binBits = binBits
    }

    output {
        File out = star_index.out
    }

}

task star_index {

    File genomeDir
    File genomeFasta
    Int threads
    Int binBits

    command {
        /usr/local/bin/STAR \
        --runThreadN ${threads} \
        --runMode genomeGenerate \
        --genomeDir ${genomeDir} \
        --genomeFastaFiles ${genomeFasta}  \
        --genomeChrBinNbits ${binBits}
    }

    runtime {
        docker: "quay.io/biocontainers/star@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
      }


    output {
        File out = genomeDir
    }

}