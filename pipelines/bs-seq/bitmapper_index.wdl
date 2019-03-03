version development

workflow bitmapper_index {
    input {
        File genome
    }

    call index{
        input: genome = genome
    }
    output {

    }
}

task index {
   input {
    File genome
   }
   command {
        /opt/BitMapperBS/bitmapperBS ~{genome}
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
  }

  output {
    File out = ""
  }
}