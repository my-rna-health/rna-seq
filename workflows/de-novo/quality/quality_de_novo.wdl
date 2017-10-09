workflow quality_de_novo {

  File reads_1
  File reads_2


  call report as initial_report_1 {
      input:
          sampleName = basename(reads_1, ".fastq.gz"),
          file = reads_1
  }

  call report as initial_report_2 {
      input:
          sampleName = basename(reads_2, ".fastq.gz"),
          file = reads_2
  }


  call trimmomatics {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        min_len = 36,
        q = 18
  }

  call report as report_trimmomatics_1 {
      input:
        sampleName = basename(trimmomatics.out1, ".fastq.gz"),
        file = trimmomatics.out1
      }

  call report as report_trimmomatics_2 {
      input:
        sampleName = basename(trimmomatics.out2, ".fastq.gz"),
        file = trimmomatics.out2
      }


  call trimmomatics_adaptive {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        min_len = 36,
        opt_len = 250
  }

  call report as report_trimmomatics_adaptive_1 {
      input:
        sampleName = basename(trimmomatics_adaptive.out1, ".fastq.gz"),
        file = trimmomatics_adaptive.out1
      }

  call report as report_trimmomatics_adaptive_2 {
      input:
        sampleName = basename(trimmomatics_adaptive.out2, ".fastq.gz"),
        file = trimmomatics_adaptive.out2
      }
}


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


task trimming_seqpurge {
  File reads_1
  File reads_2
  Int threads
  Int min_len

  command {
    SeqPurge \
      -in1 ${reads_1} \
      -in2 ${reads_2} \
      -out1 ${basename(reads_1, ".fastq.gz")}.fastq.gz \
      -out2 ${basename(reads_2, ".fastq.gz")}.fastq.gz \
      -threads ${threads} \
      -min_len ${min_len}
    }
    #-ec

    runtime {
        #docker: "quay.io/comp-bio-aging/seqpurge@sha256:5709e37c231a9fe278f23830e7f59e0806936ead10c36373256225bdf1fa4924"
        docker: "virusx/seqpurge:32d9c85395a"
    }

  output {
    File out1 = basename(reads_1, ".fastq.gz") + ".fastq.gz"
    File out2 = basename(reads_2, ".fastq.gz") + ".fastq.gz"
  }
}

task atropos_illumina_pe {
  File reads_1
  File reads_2
  Int threads

  command {
    atropos trim \
      --aligner insert \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
      -pe1 ${reads_1} \
      -pe2 ${reads_2} \
      -o ${basename(reads_1, ".fastq.gz")}_trimmed.fastq.gz \
      -p ${basename(reads_2, ".fastq.gz")}_trimmed.fastq.gz \
      --threads ${threads} \
      --correct-mismatches liberal \
      --trim-n
    }

    runtime {
        docker: "quay.io/comp-bio-aging/atropos@sha256:26624010213af4a4238a4c17c5c29e4710c178150355de945b50d1e83634854f"
    }

  output {
    File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq.gz"
    File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq.gz"
  }
}

task trimming_UrQt_pe {

  File reads_1
  File reads_2
  Int len
  Int q
  Int threads

  command {
    UrQt \
        --in ${reads_1} \
        --inpair ${reads_2} \
        --out ${basename(reads_1, ".fastq.gz")}_trimmed.fastq \
        --outpair ${basename(reads_2, ".fastq.gz")}_trimmed.fastq \
        --t ${q} --min_read_size ${len} --gz  -m ${threads}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/urqt@sha256:0c5ceb7757c5f2f6751d12861aaa08299e300012757d36c7e80551cfac3e7ba8"
  }

  output {
    File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq.gz"
    File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq.gz"
  }
}

task trimming_afterQC {

  File reads_1
  File reads_2
  Int q

  command {
    /usr/bin/pypy /opt/AfterQC/after.py \
        --read1_file ${reads_1} \
        --read2_file ${reads_2} \
        --debubble \
        -q ${q}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/afterqc@sha256:5c0bdb35ff5913f95d908c3bee0a72b0e8555597f648308e4bf6fd6391398eec"
  }

  output {
    File good = "good"
    File qc = "QC"
    File out1 = "good/" + basename(reads_1, ".fastq.gz") + ".good.fq"
    File out2 = "good/" + basename(reads_2, ".fastq.gz") + ".good.fq"
    #File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq.gz"
    #File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq.gz"
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

task trimmomatics_adaptive {

    File reads_1
    File reads_2
    Int min_len
    Int opt_len

    command {
       /usr/local/bin/trimmomatic PE \
            ${reads_1} \
            ${reads_2} \
            ${basename(reads_1, ".fastq.gz")}_trimmed.fastq.gz \
            ${basename(reads_1, ".fastq.gz")}_trimmed_unpaired.fastq.gz \
            ${basename(reads_2, ".fastq.gz")}_trimmed.fastq.gz \
            ${basename(reads_2, ".fastq.gz")}_trimmed_unpaired.fastq.gz \
            ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:1:TRUE \
            MAXINFO:${opt_len}:0.6 MINLEN:${min_len}
    }

    runtime {
        docker: "quay.io/biocontainers/trimmomatic@sha256:bf4b0b2d2747670deeb9a6adc4a50aa923b830f0b02be47e82d1b848e1368078"
    }

    output {
            File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq.gz"
            File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq.gz"
            File out1_unpaired = basename(reads_1, ".fastq.gz") + "_trimmed_unpaired.fastq.gz"
            File out2_unpaired = basename(reads_2, ".fastq.gz") + "_trimmed_unpaired.fastq.gz"
    }
}