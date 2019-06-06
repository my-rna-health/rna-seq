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

        call copy {
         input:
             files = [salmon_index.out],
             destination = indexes_folder + "/" + organism
        }
    }

    output {
        Array[File] indexes = flatten(copy.out)
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
    docker: "combinelab/salmon:0.14.0"
  }

  output {
    File out = indexName
  }

}



task copy {
    input {
        Array[File] files
        String destination
    }

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{sep=' ' files} ~{destination}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{destination}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}