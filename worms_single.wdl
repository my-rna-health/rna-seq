workflow worms_single {

  File samplesFile
  File genomeFolder

  Array[Array[File]] samples = read_tsv(samplesFile)

  scatter (sample in samples) {
    call extract_single_fastq {
        input:
            fileName = sample[0],
            file =  sample[1]
    }

    call report as initial_report {
        input:
            fileName = sample[0],
            file = extract_single_fastq.out
    }

    call adapter_removal_single{
        input:
            fileName = sample[0],
            file = extract_single_fastq.out
    }

    call report as trimming_report {
            input:
                fileName = sample[0],
                file =  adapter_removal_single.out
        }

    call star {
                input:
                    file =  adapter_removal_single.out,
                    genomeDir = genomeFolder
            }
  }
}


task stats {

  String fileName
  File file


  command {
    /opt/sratoolkit/sra-stat ${file} > ${fileName}_stats.txt
  }

  runtime {
    docker: "itsjeffreyy/sratoolkit@sha256:9938a78b61b702992e28202e60c60e84ede9d6495b84edd551f6c3e9d504d33d"
  }

  output {
    File stats = "${fileName}_stats.txt"
  }
}

task extract_single_fastq {

  String fileName
  File file

  # read the following explanations for parameters
  # https://edwards.sdsu.edu/research/fastq-dump/

  #command {
  #  fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-files --clip ${file}
  #}

  command {
    /opt/sratoolkit/fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-files --clip ${file}
  }

  runtime {
    #docker: "quay.io/biocontainers/sra-tools:2.8.1--0"
    docker: "itsjeffreyy/sratoolkit@sha256:9938a78b61b702992e28202e60c60e84ede9d6495b84edd551f6c3e9d504d33d"
  }

  output {
    File out = "${fileName}_pass_1.fastq.gz"
  }

}

task extract_sam {

  String fileName
  File file

  command {
    sam-dump ${file} > ${fileName}.sam
  }

  runtime {
    #docker: "quay.io/biocontainers/sra-tools:2.8.1--0"
    #docker: "itsjeffreyy/sratoolkit"
    docker: "itsjeffreyy/sratoolkit@sha256:9938a78b61b702992e28202e60c60e84ede9d6495b84edd551f6c3e9d504d33d"
  }

  output {
    File sam = "${fileName}.sam"
  }
}

task report {

  String fileName
  File file

  command {
    /opt/FastQC/fastqc ${file}
  }

  runtime {
    docker: "quay.io/ucsc_cgl/fastqc@sha256:86d82e95a8e1bff48d95daf94ad1190d9c38283c8c5ad848b4a498f19ca94bfa"
    #docker: "quay.io/biocontainers/fastqc@sha256:bb57a4deeec90633e746afbc38c36fdb202599fe71f9557b94652e9c8f3c1a02"
  }

  output {
    #File out = sub(file, "\\.fastq.gz", "_fastqc.gz")
    File out = "${fileName}.fastqc.gz"
  }
}

task adapter_removal_single {

  String fileName
  File file

  command {
    sickle se -g -q 22 -f ${file} -t sanger -o ${fileName}"_trimmed.fastq.gz"
  }

  runtime {
    #docker: "ksimons77/sickle@sha256:3e3fed0fe7735e47998c1b82b8bb920a542530479041c48aa2fe21ca8f9ee0a3"
    docker: "asidorovj/sickle@sha256:933e4a880c58804248179c3819bb179c45ba86c85086d8435d8ab6cf82bca63c"
  }

  output {
    File out = "${fileName}_trimmed.fastq.gz"
  }
}

task cutadapt {
String fileName
  File file

  command {
    sickle se -g -q 22 -f ${file} -t sanger -o ${fileName}"_trimmed.fastq.gz"
  }

  runtime {
    docker: "jdidion/atropos@sha256:dbdf89b1fb65c34d2a7fc3483030599b08556ba96a37b888b6f1595c023cf995"
  }

  output {
    File out = "${fileName}_trimmed.fastq.gz"
  }

}



task star {

  Int  numberOfThreads = 8
  File file
  File genomeDir

  command {
    STAR --runThreadN ${numberOfThreads} --genomeDir ${genomeDir} --readFilesCommand gunzip -c --readFilesIn ${file}
  }

  runtime {
    docker: "quay.io/biocontainers/star@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
  }

  output {
    String result = "STAR WORKS!"
  }

}
