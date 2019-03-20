workflow quant_index {

  Map[String, File] transcriptomes
  String results_folder = "/data/indexes/salmon"

   scatter(kv in transcriptomes){
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

}

task salmon_index {

  File transcriptomeFile
  String indexName

  command {
    salmon --no-version-check index -t ${transcriptomeFile} -i ${indexName} --type quasi
  }

  runtime {
    docker: "combinelab/salmon:0.12.0"
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