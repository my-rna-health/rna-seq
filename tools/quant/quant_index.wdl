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
    salmon --no-version-check index -t ${transcriptomeFile} -i ${indexName} --type quasi
  }

  runtime {
    docker: "combinelab/salmon@sha256:bb9b64804d9ac79c98cc19c11a61e65bb290446beec377d46229c2686990c311" #0.11.2
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