workflow quant_index {

  File transcriptome
  String name
  String indexes_folder #will be created if needed

  call salmon_index  {
      input:
          transcriptomeFile = transcriptome,
          indexName =  name
  }

  call copy {
    input:
        files = [salmon_index.out],
        destination = indexes_folder
  }

  output {
    Array[File] out = copy.out
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