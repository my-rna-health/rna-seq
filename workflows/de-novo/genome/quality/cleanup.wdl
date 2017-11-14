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


  #quay.io/comp-bio-aging/geoparse  --location /data --filetype sra --keep_sra true GSM1696283 GSM1696284
  command {
    /opt/geoparse/run.py --location ./ --filetype fastq --keep_sra false ${sample}
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