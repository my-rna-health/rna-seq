workflow assembly_trinity {


  String max_memory
  Int threads
  File alignment

  call trinity_assembly {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        cores = threads,
        max_memory = max_memory
  }

}

task trinity_guided_assembly {

    File aligment_bam
    Int cores
    String max_memory
    File alignment

    command {
        Trinity --genome_guided_bam ${aligment_bam} \
                  --max_memory ${max_memory} \
                  --CPU ${cores}

    }

    runtime {
        docker: "trinityrnaseq/trinityrnaseq@sha256:9f607b309e9c2d3457775ad21a54179e51cafbc76d4f1af97ca8ef218b128637"
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
}