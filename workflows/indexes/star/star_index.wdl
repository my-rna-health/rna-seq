workflow StarIndex {
    File referenceGenome
    File indexDir
    Int threads
    Int binBits

    call star_index {
        input:
            genomeDir = indexDir,
            genomeFasta = referenceGenome,
            threads = threads,
            binBits = binBits
    }

    call copy_dir {
        input:
            source = star_index.out,
            destination = indexDir
    }


    output {
        File out = copy_dir.out
    }

}

task star_index {

    File genomeDir
    File genomeFasta
    Int threads
    Int binBits
    
    command {
        mkdir INDEX
        /usr/local/bin/STAR \
        --runThreadN ${threads} \
        --runMode genomeGenerate \
        --genomeDir INDEX \
        --genomeFastaFiles ${genomeFasta}  \
        --genomeChrBinNbits ${binBits} \
        --limitGenomeGenerateRAM=100000000000
    }

    runtime {
        docker: "quay.io/biocontainers/star@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
      }


    output {
        File out = "INDEX"
    }

}

task copy_dir {
    File source
    File destination

    command {
        mkdir -p ${destination}
        cp -R -H -u ${source}/* ${destination}
    }

    output {
        File out = destination
    }
}

task copy {
    Array[File] files
    File destination

    command {
        mkdir -p ${destination}
        cp -R -H -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}