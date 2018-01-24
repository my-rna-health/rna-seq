workflow download {

    String sample
    String filetype
    Boolean keep_sra

    #call get_sample{
    #input:
    #    sample = sample,
    #    filetype = filetype,
    #    keep_sra = keep_sra
    #}

    call fastq_dump{
        input:
            sample = sample
    }
}

task fastq_dump {
    String sample

    # read the following explanations for parameters
    # https://edwards.sdsu.edu/research/fastq-dump/

    command {
        fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-files --clip ${sample}
    }

    runtime {
        docker: "quay.io/biocontainers/sra-tools@sha256:6d8c1daefedf6a4d00b3351d030228ca7cc4669f73f5fb76d76d1b80e79905f1"
    }

    output {
        File reads_1 = sample + "_1.fastq"
        File reads_1 = sample + "_2.fastq"
    }
}

task get_sample {

    String sample
    String filetype
    Boolean keep_sra

    # read the following explanations for parameters
    # https://edwards.sdsu.edu/research/fastq-dump/

    #command {
    #  fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-files --clip ${file}
    #}

  #quay.io/comp-bio-aging/geoparse  --location /data --filetype sra --keep_sra true GSM1696283 GSM1696284
  command {
    /opt/geoparse/run.py --location ./ --filetype ${filetype} --keep_sra ${keep_sra} ${sample}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/geoparse@sha256:99cdeeaac6b3e3060a8d088e99988d9b19ba76ae6aec8f2e3e71218313361f2a"
  }

  output {
    Array[Array[String]] out = read_tsv("output.tsv")
    String sampleName = basename(out[0][0])
    File outputFile = out[0][1]
  }

}

