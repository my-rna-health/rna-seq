workflow quantify {

  File transcriptome
  Array[File] samples

  scatter (sample in samples) {
      call salmon as quantification  {
          input:
              index = transcriptome,
              file =  sample
      }
  }

  output {
    Array[File] out = quantification.out
  }

}

task salmon {
  File file
  File index

  command {
    salmon quant -i ${index} -l A -r ${file} -o transcripts_quant
  }

  runtime {
    docker: "combinelab/salmon@sha256:2d09e1113f5bf1aa6be9354b9c8be556a059bbb3cdc46067080894b6ebb2a983"
  }

  output {
    File out = "transcripts_quant"
  }
}

task kallisto {

  File index
  File file
  Int len

  command {
    kallisto quant -i ${index} -o output --single -l ${len} ${file}
  }

  runtime {
    docker: "quay.io/ucsc_cgl/kallisto@sha256:83712616f897deccae1266050eb309835a177be95dfa2c8de6381ed16c2f55ab"
  }

  output {
    File out = "output"
  }

}

task sleuth {

    command {

    }

    runtime {
        docker: "quay.io/biocontainers/r-sleuth"
    }

    output {
        String out = "it works!"
    }

}
