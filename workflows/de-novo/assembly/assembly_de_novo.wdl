workflow assembly_de_novo {

  File reads_1
  File reads_2

  call trimmomatics {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        min_len = 36,
        q = 19
  }

  call trinity_assembly {
      input:
        reads_1 = trimmomatics.out1,
        reads_2 = trimmomatics.out2,
        cores = 2,
        max_memory = "32G"
  }

}

task trimmomatics {
    File reads_1
    File reads_2
    Int q
    Int min_len

    command {
       /usr/local/bin/trimmomatic PE \
            ${reads_1} \
            ${reads_2} \
            ${basename(reads_1, ".fastq.gz")}_trimmed.fastq.gz \
            ${basename(reads_1, ".fastq.gz")}_trimmed_unpaired.fastq.gz \
            ${basename(reads_2, ".fastq.gz")}_trimmed.fastq.gz \
            ${basename(reads_2, ".fastq.gz")}_trimmed_unpaired.fastq.gz \
            ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:1:TRUE \
            SLIDINGWINDOW:4:${q} MINLEN:${min_len}
    }

    runtime {
        docker: "quay.io/biocontainers/trimmomatic@sha256:bf4b0b2d2747670deeb9a6adc4a50aa923b830f0b02be47e82d1b848e1368078"
    }

    output {
        File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq.gz"
        File out1_unpaired = basename(reads_1, ".fastq.gz") + "_trimmed_unpaired.fastq.gz"
        File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq.gz"
        File out2_unpaired = basename(reads_2, ".fastq.gz") + "_trimmed_unpaired.fastq.gz"

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
        File out_dir = "trinity_out_dir"
        File out = "trinity_out_dir/Trinity.fasta"
    }
}
