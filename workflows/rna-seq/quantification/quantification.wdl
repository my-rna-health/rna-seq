workflow quantification {

    File batch
    File references
    Int threads
    File samples_folder
    String results_folder

    Array[Array[String]] samples = read_tsv(batch)
    #Map[String, Map[String, String]] indexes = read_json(references)
    #Map[String, Map[String, String]] indexes = read_json(references)

    #Object indexes = read_json(references)


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