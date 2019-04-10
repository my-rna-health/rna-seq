version development

workflow quant_index {
    input {
        Map[String, File] transcriptomes
        String results_folder = "/data/indexes/salmon"
    }
    Array[Pair[String, File]] ts = transcriptomes


    scatter (kv in ts) {
        String organism = sub(kv.left, " ", "_")
        File transcriptome = kv.right

        call salmon_index  {
           input:
               transcriptomeFile = transcriptome,
               indexName =  organism
        }
    }
    call copy {
     input:
         files = salmon_index.out,
         destination = results_folder
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
    docker: "combinelab/salmon:0.13.1"
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