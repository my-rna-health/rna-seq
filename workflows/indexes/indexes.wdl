workflow indexes {

    File indexesFolder
    String genome #.fa.gz
    String transcriptome #.fa.gz
    String annotation #.gtf
    String species #species and also name of the index/
    String release #release version
    String currentIndexFolder = indexesFolder + "/" + species + "/" + release

    call download as download_genome {
        input:
             url = genome,
             folder = currentIndexFolder,
             filename = basename(genome)
    }

    call download as download_transcriptome {
        input:
             url = transcriptome,
             folder = currentIndexFolder,
             filename = basename(transcriptome)
    }

    call download as download_annotation {
        input:
             url = annotation,
             folder = currentIndexFolder,
             filename = basename(annotation)
    }

    call salmon_index {
        input:
            transcriptomeFile = download_transcriptome.out,
            indexName = species
    }

    call copy {
        input:
            file = salmon_index.out,
            folder = currentIndexFolder + "/" + "salmon"
    }

    output {
        File out = copy.out
    }


}

task download {
    String url
    String folder
    String filename

    command {
        mkdir -p ${folder}
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 -nc -q -O "${folder}/${filename}" "${url}"
      }

    output {
        File out = folder + "/" + filename
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
    File out = indexName
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

task move {
    File file
    String folder # String as the folder may not exist

    command {
        mkdir -p ${folder} && mv -u ${file} ${folder}
    }

  output {
    File out = folder + "/" + basename(file)
  }
}

task copy {
    File file
    String folder # String as the folder may not exist

    command {
        mkdir -p ${folder} && cp -R -u ${file} ${folder}
    }

  output {
    File out = folder + "/" + basename(file)
  }
}