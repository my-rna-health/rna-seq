workflow assembly_trinity {


  String max_memory
  Int threads
  File alignment

  call trinity_guided_assembly {
      input:
        aligment_bam = alignment,
        cores = threads,
        max_memory = max_memory
  }

}

task trinity_guided_assembly {

    File aligment_bam
    Int cores
    String max_memory
    Int max_intron = 350000

    command {
        Trinity --genome_guided_bam ${aligment_bam} \
                  --max_memory ${max_memory} \
                  --genome_guided_max_intron ${max_intron} \
                  --bflyCalculateCPU \
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