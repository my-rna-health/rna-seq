task report {

  String sampleName
  File file

  command {
    /opt/FastQC/fastqc ${file} -o .
  }

  runtime {
    docker: "quay.io/ucsc_cgl/fastqc@sha256:86d82e95a8e1bff48d95daf94ad1190d9c38283c8c5ad848b4a498f19ca94bfa"
  }

  output {
    File out = sampleName+"_fastqc.zip"
  }
}


task trimming_sickle_pe {

  File reads_1
  File reads_2
  Int len
  Int q

  command {
    sickle pe \
            -f ${reads_1} \
            -r ${reads_2} \
            -t sanger \
            -o ${basename(reads_1, ".fastq.gz")}_trimmed.fastq \
            -p ${basename(reads_2, ".fastq.gz")}_trimmed.fastq \
            -s unpaired.fastq \
            -q ${q} -l ${len}
  }

  runtime {
    docker: "asidorovj/sickle@sha256:933e4a880c58804248179c3819bb179c45ba86c85086d8435d8ab6cf82bca63c"
  }

  output {
    File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq"
    File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq"
  }
}

task trimming_skewer {

  File reads_1
  File reads_2
  Int len
  Int q

  command {
    skewer pe \
            -m pe -q ${q}  -n -u -l ${len} \
            ${reads_1} ${reads_2}
  }

  runtime {
    docker: "quay.io/biocontainers/skewer@sha256:047a72bb4dc61d9896318beb67f90e71bf2557c54bdd1142cea8820e516607a1"
  }

  output {
    File out1 = basename(reads_1, ".fastq.gz") + ".trimmed.fastq"
    File out2 = basename(reads_2, ".fastq.gz") + ".trimmed.fastq"
  }
}

task trimming_adapter_removal {

  File reads_1
  File reads_2
  Int q
  Int threads

  command {
    AdapterRemoval \
        --file1 ${reads_1} \
        --file2 ${reads_2} \
        --output1 ${basename(reads_1, ".fastq.gz")}_trimmed.fastq \
        --output2 ${basename(reads_2, ".fastq.gz")}_trimmed.fastq \
        --threads ${threads} --trimwindows ${q} \
        --trimns --trimqualities --collapse
  }

  runtime {
    docker: "quay.io/biocontainers/skewer"
  }

  output {
    File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq"
    File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq"
  }
}
