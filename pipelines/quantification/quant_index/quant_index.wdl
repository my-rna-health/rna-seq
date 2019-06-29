version development

struct Transcriptome{
    String species
    String version
    File reference
}

workflow quant_index {
    input {
        Array[Transcriptome] transcriptomes
        String indexes_folder
    }

    scatter (trans in transcriptomes) {
        String organism = sub(trans.species, " ", "_")
        call salmon_index  {
           input:
               transcriptomeFile = trans.reference,
               indexName =  trans.version
        }

        call copy_folder {
         input:
             dir = salmon_index.out,
             destination = indexes_folder + "/" + organism
        }
    }

    output {
        Array[File] indexes = copy_folder.out
    }

}

task salmon_index {
    input {
        File transcriptomeFile
        String indexName
    }

  command {
    salmon --no-version-check index -t ~{transcriptomeFile} -i ~{indexName} --type quasi
  }

  runtime {
    docker: "combinelab/salmon:0.14.1"
  }

  output {
    Directory out = indexName
  }

}

task copy_folder {
    input {
        Directory dir
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{dir} ~{where}
    }

    output {
        File out = where
    }
}