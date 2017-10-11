workflow assembly_trinity {

  File reads_1
  File reads_2

  call trinity_assembly {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        cores = 2,
        max_memory = "1G"
  }

}

task trinity_assembly {

    File reads_1
    File reads_2
    Int cores
    String max_memory

    command {
        Trinity --seqType fq  \
         --left ${reads_1} \
         --right ${reads_2}  \
         --max_memory ${max_memory} \
         --CPU ${cores}
    }

    runtime {
        docker: "trinityrnaseq/trinityrnaseq@sha256:9f607b309e9c2d3457775ad21a54179e51cafbc76d4f1af97ca8ef218b128637"
    }

    output {


    }
}