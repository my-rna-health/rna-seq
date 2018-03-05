workflow quantification {

    File batch
    File references
    Int threads
    File samples_folder
    String results_folder

    call prepare_samples {
        input: samples = batch,references = references, samples_folder = samples_folder
    }

    call copy as report_invalid {
        input: files = [prepare_samples.invalid], destination = results_folder + "/invalid"
    }


    scatter(row in read_tsv(prepare_samples.novel)) {
            #"GSM",	"GSE",	"Species",	"Sequencer",
            #"Type", "Sex",	"Age",	"Tissue",
            #"Extracted molecule", "Strain",
            #"Comments", "salmon", "transcriptome", "gtf"
        String gsm = row[0]
        String gse = row[1]

        call echo {
            input: message = gsm
        }

    }

}

task echo {
    String message

    command {
        echo ${message} >> /pipelines/test/echo.txt
    }

    output {
        String out = message
    }
}


task prepare_samples {
    File samples
    File references
    File samples_folder

    command {
        /scripts/run.sc --samples ${samples} --references ${references} --cache ${samples_folder}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/prepare-samples:latest"
    }

    output {
        File invalid = "invalid.tsv"
        File novel = "novel.tsv"
        File cached = "cached.tsv"
    }
}


task print {
    String incoming
    String where

    command {
        echo ${incoming} > ${where}
    }

    output {
        String out =  incoming
    }
}

task get_sample {

  String sample

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

task salmon {
  File index
  File reads_1
  File reads_2
  Int numThreads

  command {
    salmon quant -i ${index} --threads ${numThreads} -l A -1 ${reads_1} -2 ${reads_2} -o transcripts_quant
  }

  runtime {
    docker: "combinelab/salmon@sha256:8a5f0de02b0df1b2571f8200e276c09ef1dd499ca13a883577230d85d8e644c3"
  }

  output {
    File out = "transcripts_quant"
  }
}


task copy {
    Array[File] files
    String destination

    command {
        mkdir -p ${destination}
        cp -L -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}