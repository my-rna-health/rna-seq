version development
workflow salmon_merge {
    input {
        String destination
        Array[File] transcripts
    }
}

task merge_expressions {
    input{
        Int p = 3
        String max_memory = "20G"
    }


  command {
    salmon index -t ~{gentrome} -d ~{decoys} -p ~{p} -i ~{indexName}
  }

  runtime {
    docker: "quay.io/biocontainers/salmon@sha256:e56485bfa26913aebaa6351b2ddb1308d0dc0352bf15e7f5431bc58ba5465809" #1.9.0--h7e5ed60_1
    maxRetries: 3
    docker_memory: "${max_memory}"
    docker_cpu: "${p}"
  }

  output {

  }
}

task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{where}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}