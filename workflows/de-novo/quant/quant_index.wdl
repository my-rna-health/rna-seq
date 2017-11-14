workflow quant_index {

  File transcriptome
  String name

  call salmon_index  {
      input:
          transcriptomeFile = transcriptome,
          indexName =  name
  }

  output {
    File out = salmon_index.out
  }

}

task salmon_index {

  File transcriptomeFile
  String indexName

  command {
    salmon index -t ${transcriptomeFile} -i ${indexName} --type quasi
  }

  runtime {
    docker: "combinelab/salmon@sha256:2d09e1113f5bf1aa6be9354b9c8be556a059bbb3cdc46067080894b6ebb2a983"
  }

  output {
    File out = indexName
  }

}