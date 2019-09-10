version development

workflow bitmapper_index {

    input {
        File genome
        String index_name
        String destination
    }
    call index{
                input: genome = genome,
                index_folder = index_name
    }
    call copy {
        input:
            files = [index.out],
            destination = destination
    }

    output {
        Array[File] out = copy.out
    }
}

task index {
   input {
        File genome
        String index_folder
   }

   command {
        /opt/BitMapperBS/bitmapperBS --index ~{genome}  --index_folder ~{index_folder}
   }

   output {
        File out = index_folder
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
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