workflow tester {

  File samplesFile
  File genomeFolder

  Array[Array[File]] samples = read_tsv(samplesFile)

  scatter (sample in samples) {

    call extract_single_fastq {
        input:
            sampleName = sample[0],
            file =  sample[1]
    }

    call detect_adapters {
        input:
            sampleName = sample[0],
            file = extract_single_fastq.out
    }


    call adapter_trimming_atropos{
        input:
            sampleName = sample[0],
            file = extract_single_fastq.out,
            adapters = detect_adapters.out
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

task detect_adapters {
  String sampleName
  File file

  command {
    /usr/local/bin/atropos detect -se ${file} -d heuristic -O fasta -o ${sampleName}_adapters.fasta
  }

  runtime {
    docker: "quay.io/comp-bio-aging/atropos@sha256:e83e9c13107649c5d56a6dd3d5573919cd7a9a4b7048b894c503fd44d4f29a68"
  }

  output {
    File out = "${sampleName}_adapters.0.fasta"
  }

}

task adapter_trimming_atropos {

  String sampleName
  File file
  File adapters

  command {
    /usr/local/bin/atropos trim --trim-n -se ${file} --known-adapters-file ${adapters} -o ${sampleName}"_trimmed.fastq.gz"
  }

  runtime {
    docker: "quay.io/comp-bio-aging/atropos@sha256:e83e9c13107649c5d56a6dd3d5573919cd7a9a4b7048b894c503fd44d4f29a68"
  }

  output {
    File out = "${sampleName}_trimmed.fastq.gz"
  }

}