workflow download {

    String sample
    String filetype
    String keep_sra

    call get_sample{
    input:
        sample = sample,
        filetype = filetype,
        keep_sra = keep_sra
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
    /opt/geoparsec/run.py --location ./ --filetype ${filetype} --keep_sra false ${sample}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/geoparse@sha256:c019939664836a95c6e715d0d47dac479963ddfee821e28452dff3fcd5e41b1b"
  }

  output {
    Array[Array[String]] out = read_tsv("output.tsv")
    String sampleName = basename(out[0][0])
    File outputFile = out[0][1]
  }

}
