workflow assembly_de_novo {

  File reads_1
  File reads_2

  call trinity_assembly {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        cores = 8
  }

}

task trinity_assembly {

    File reads_1
    File reads_2
    Int cores


#TODO: fix
    command {
        Trinity --seqType fq --max_memory 50G  \
         --left ${reads_1} \
         --right ${reads_1}  \
         --CPU ${cores}
    }

    runtime {
        docker: "trinityrnaseq/trinityrnaseq@sha256:9f607b309e9c2d3457775ad21a54179e51cafbc76d4f1af97ca8ef218b128637"
    }

}