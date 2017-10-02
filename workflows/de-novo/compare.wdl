workflow compare {

   File reads_1
   File reads_2

}

task report {

  String sampleName
  File file

  command {
    /opt/FastQC/fastqc ${file} -o .
  }

  runtime {
    docker: "quay.io/ucsc_cgl/fastqc@sha256:86d82e95a8e1bff48d95daf94ad1190d9c38283c8c5ad848b4a498f19ca94bfa"
    #docker: "quay.io/biocontainers/fastqc@sha256:bb57a4deeec90633e746afbc38c36fdb202599fe71f9557b94652e9c8f3c1a02"
  }

  output {
    #File out = sub(file, "\\.fastq.gz", "_fastqc.gz")
    File out = "${sampleName}_fastqc.zip"
  }
}

task trimmomatic {

   File reads_1
   File reads_2

    command {
        trimmomatic ${reads_1} ${reads_2} SLIDINGWINDOW:4:20 MINLEN:36
    }

    runtime {
        docker: "quay.io/biocontainers/trimmomatic@sha256:bf4b0b2d2747670deeb9a6adc4a50aa923b830f0b02be47e82d1b848e1368078"
    }
}

task adapter_removal {

  File reads_1
  File reads_2

  command {
    AdapterRemoval --file1 ${reads_1} --file2 ${reads_2} --basename output_paired --trimns --trimqualities --collapse
  }

   runtime {
      docker: "quay.io/ucsc_cgl/fastqc@sha256:86d82e95a8e1bff48d95daf94ad1190d9c38283c8c5ad848b4a498f19ca94bfa"
      #docker: "quay.io/biocontainers/fastqc@sha256:bb57a4deeec90633e746afbc38c36fdb202599fe71f9557b94652e9c8f3c1a02"
    }
}

task adapter_trimming_sickle {

  String sampleName
  File file

  command {
    /usr/local/bin/sickle/sickle se -f ${file} -t sanger  -q 25 -g -o "${sampleName}_trimmed.fastq.gz"
  }

  runtime {
    #docker: "ksimons77/sickle@sha256:3e3fed0fe7735e47998c1b82b8bb920a542530479041c48aa2fe21ca8f9ee0a3"
    docker: "asidorovj/sickle@sha256:933e4a880c58804248179c3819bb179c45ba86c85086d8435d8ab6cf82bca63c"
  }

  output {
    File out = "${sampleName}_trimmed.fastq.gz"
  }
}

