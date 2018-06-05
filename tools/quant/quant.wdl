workflow quantification {

    File index
    File reads_1
    File reads_2
    Int threads
    String results_folder


    call salmon {
        input:
            index = index,
            numThreads = threads,
            reads_1 = reads_1,
            reads_2 = reads_2
    }

    call copy {
        input:
            destination = results_folder,
            files = [salmon.out]
    }

}

task salmon {
  File index
  File reads_1
  File reads_2
  Int numThreads

  command {
    salmon quant -i ${index} --threads ${numThreads} -l A -1 ${reads_1} -2 ${reads_2} -o transcripts_quant
  }

  runtime {
    docker: "combinelab/salmon@sha256:031d53d3da93887acec49a4fd0a4d4776cb9057acefa6fedf3faf655ab7bab4a"
  }

  output {
    File out = "transcripts_quant"
  }
}


task copy {
    Array[File] files
    String destination

    command {
        mkdir -p ${destination}
        cp -L -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}