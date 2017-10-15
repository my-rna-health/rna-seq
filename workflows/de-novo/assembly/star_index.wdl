workflow StarIndex {
    File referenceGenome
    File indexDir
    Int threads
    Int binBits
    String resultsFolder #will be created if needed

    call star_index {
        input:
            genomeDir = indexDir,
            genomeFasta = referenceGenome,
            threads = threads,
            binBits = binBits
    }

    call copy as copy_results {
        input:
            files = [star_index.out],
            destination = results_folder
    }


    output {
        File out = copy_results.out
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
        --genomeChrBinNbits ${binBits} \
        --limitGenomeGenerateRAM=100000000000
    }

    runtime {
        docker: "quay.io/biocontainers/star@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
      }


    output {
        File out = genomeDir
    }

}

task copy {
    Array[File] files
    File destination

    command {
        mkdir -p ${destination}
        cp -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}