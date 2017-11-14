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
