version development

workflow bs_map {
    input {
        File genome_index
        Array[File] reads
        String filename
        String destination
        Boolean is_paired = true
        Int map_threads = 8
    }


    call bitmapper_map {
        input:
            index_folder = genome_index,
            reads = reads,
            is_paired = is_paired,
            threads = map_threads
    }

    call bitmapper_extract {
        input:
            index_folder = genome_index,
            bmms = bitmapper_map.out,
            filename = filename
    }

    call copy as copy_bit {
        input:
            files = [bitmapper_extract.out],
            destination = destination
    }

}

task bitmapper_map {
   input {
        File index_folder
        Array[File] reads
        Boolean is_paired
        Int threads
   }
   command {
        /opt/BitMapperBS/bitmapperBS --search ~{index_folder} ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]} -t ~{threads} --methy_out
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
  }

  output {
    Array[File] out = glob("*.bmm")
  }
}

task bitmapper_extract {
   input {
        File index_folder
        Array[File] bmms
        String filename
   }
   command {
        /opt/BitMapperBS/bitmapperBS --methy_extract ~{index_folder} --seq ~{bmms} --output ~{filename}.bedGraph
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
  }

  output {
    File out = filename + ".bedGraph"
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