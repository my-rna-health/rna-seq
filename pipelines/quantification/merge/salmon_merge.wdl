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
    docker: "quay.io/biocontainers/salmon@sha256:0aea466ba3eae62cb3ea8077e30b6212ca152f73e9520ca19f7421c5be519ef9" #1.1.0--hf69c8f4_0
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