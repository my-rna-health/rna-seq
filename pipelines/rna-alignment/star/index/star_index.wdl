version development

workflow StarIndex {

    input {
        File referenceGenome
        File indexDir
        File gtf
        Int threads = 16
    }

     call make_folder{
        input: path = indexDir
     }

    call star_index {
        input:
            genomeDir = make_folder.out,
            genomeFasta = referenceGenome,
            threads = threads,
            gtf = gtf
    }


    output {
        File out = star_index.out
    }
}

task make_folder{
    input {
        String path
    }

    command {
        mkdir -p ~{path}
    }

    output {
        File out = path
    }
}

task star_index {
    input {
        File genomeDir
        File genomeFasta
        File gtf
        Int threads
    }

    command {
        /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --runMode genomeGenerate \
        --genomeDir ~{genomeDir} \
        --genomeFastaFiles ~{genomeFasta}  \
        --sjdbGTFfile ~{gtf}
        --limitGenomeGenerateRAM=58000000000
    }

    runtime {
        docker: "quay.io/biocontainers/star@sha256:6556f4d3a9f767f93cd8841ef0ac0f15d29d2c69b6a6f1f9b6ccc3e1da40207b" #2.7.1a--0
      }


    output {
        File out = genomeDir
    }

}

task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{where}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}