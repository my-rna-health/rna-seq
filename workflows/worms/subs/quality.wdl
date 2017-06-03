workflow cleanup {

  Array[Array[File]] samples

  scatter (sample in samples) {

    call extract_single_fastq {
        input:
            sampleName = sample[0],
            file =  sample[1]
    }

    call report as initial_report {
        input:
            sampleName = "${sample[0]}_pass_1",
            file = extract_single_fastq.out
    }

    call adapter_trimming_sickle {
        input:
            sampleName = sample[0],
            file = extract_single_fastq.out
    }

    call report as trimming_report {
            input:
                sampleName = "${sample[0]}_trimmed",
                file = adapter_trimming_sickle.out
        }
  }

  output {
        Map[String, Array[File]] out = {
            "trimmed_report": trimming_report.out,
            "trimmed": adapter_trimming_sickle.out,
            "initial_report": initial_report.out
        }
  }

}

task extract_single_fastq {

  String sampleName
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
    File out = "${sampleName}_pass_1.fastq.gz"
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
    #docker: "quay.io/biocontainers/fastqc@sha256:bb57a4deeec90633e746afbc38c36fdb202599fe71f9557b94652e9c8f3c1a02"
  }

  output {
    #File out = sub(file, "\\.fastq.gz", "_fastqc.gz")
    File out = "${sampleName}_fastqc.zip"
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

