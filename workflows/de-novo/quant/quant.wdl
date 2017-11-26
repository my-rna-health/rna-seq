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
    docker: "combinelab/salmon@sha256:2d09e1113f5bf1aa6be9354b9c8be556a059bbb3cdc46067080894b6ebb2a983"
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