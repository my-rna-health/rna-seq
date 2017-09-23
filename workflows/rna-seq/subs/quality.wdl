workflow cleanup {

    Array[String] condition

    String condition_name = condition[0]

    scatter (sample in condition) {
        #first sample is not a sample but name of the condition
        if(sample != condition_name) {
            call get_sample {input: condition = condition_name, sample = sample }
        }
    }

    scatter (sample in select_all(get_sample.out)) {
        call report as initial_report {
            input:
                sampleName = get_sample.sampleName,
                file = get_sample.outputFile
        }

        call adapter_trimming_sickle {
            input:
                sampleName = get_sample.sampleName,
                file = get_sample.outputFile
        }

        call report as trimming_report {
            input:
                sampleName = basename(adapter_trimming_sickle.out),
                file = adapter_trimming_sickle.out
                }
    }

  output {
        Map[String, Array[File]] out = {
            "trimmed_report": select_all(trimming_report.out),
            "trimmed": select_all(adapter_trimming_sickle.out),
            "initial_report": select_all(initial_report.out)
        }
  }

}

task get_sample {

    String condition
    String sample

    # read the following explanations for parameters
    # https://edwards.sdsu.edu/research/fastq-dump/

    #command {
    #  fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-files --clip ${file}
    #}

  #quay.io/comp-bio-aging/geoparse  --location /data --filetype sra --keep_sra true GSM1696283 GSM1696284
  command {
    /opt/geoparsec/run.py --location ./ --filetype fastq --keep_sra false ${sample}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/geoparse@sha256:174e51df9dc32166b0675afc2f3d4a73d44f6a69c3cfae9ed7d52367a9cc4222"
  }

  output {
    Array[Array[String]] out = read_tsv("output.tsv")
    String sampleName = basename(out[0][0])
    File outputFile = out[0][1]
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

