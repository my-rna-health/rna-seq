workflow indexes {

  File genomeFolder
  File transcriptomeFile
  String indexName

 call salmon_index {
    input:
        transcriptomeFile = transcriptomeFile,
        indexName = indexName
 }

  output {
    File out = salmon_index.out
  }

}

task kallisto {

  File file
  String indexName

  command {
    kallisto index ${file} -i ${indexName}
  }

  runtime {
    docker: "quay.io/ucsc_cgl/kallisto@sha256:83712616f897deccae1266050eb309835a177be95dfa2c8de6381ed16c2f55ab"
  }

  output {
    File out = "${indexName}"
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
    File out = "${indexName}"
  }

}