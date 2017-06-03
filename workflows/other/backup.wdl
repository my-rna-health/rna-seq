
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
            file = extract_single_fastq.out
    }

    call adapter_removal_single{
        input:
            file = extract_single_fastq.out
    }

    call report as trimming_report {
            input:
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
    /opt/sratoolkit/sra-stat ${file} > stats.txt
  }

  runtime {
    docker: "itsjeffreyy/sratoolkit"
  }

  output {
    File stats = "stats.txt"
  }
}

task extract {

  String fileName
  File file

  command {
    /opt/sratoolkit/sam-dump ${file} > stats.txt
  }

  runtime {
    docker: "itsjeffreyy/sratoolkit"
  }

  output {
    File stats = "${fileName}.sam"
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



task adapter_trimming_sickle {

  File file

  File result =  sub(file, "\\.fastq$", "_trimmed.fastqc.gz")

  command {
    /usr/local/bin/sickle/sickle se -f ${file} -t sanger  -q 25 -g ${result}
  }

  runtime {
    #docker: "ksimons77/sickle@sha256:3e3fed0fe7735e47998c1b82b8bb920a542530479041c48aa2fe21ca8f9ee0a3"
    docker: "asidorovj/sickle@sha256:933e4a880c58804248179c3819bb179c45ba86c85086d8435d8ab6cf82bca63c"
  }

  output {
    File out = sub(file, "\\.fastq$", "_trimmed.fastqc.gz")
  }
}

task adapter_trimming_trim_galore {
  String fileName
  File file

  #https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

  command {
    trim_galore --fastqc --illumina --trim-n ${file}
  }

  runtime {
    docker: "quay.io/biocontainers/trim-galore@sha256:acba2c901dc982d79615ac1d71f2c2ab444b844d976d6308f4d7c9ca06bdb025"
  }

  output {
    File out = "${fileName}_trimmed.fastq.gz"
  }
}


task star_index {

  #TODO: finish STAR part
  Int  numberOfThreads = 8
  File genomeDir
  File genomeFasta
  File annotation

  command {
    STAR --runThreadN ${numberOfThreads} --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles ${genomeFasta} --sjdbGTFfile ${annotation}
  }

  runtime {
    docker: "quay.io/biocontainers/star:2.5.3a--0@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
  }

  output {
    File out = "output_single.fastq.gz"
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
    docker: "quay.io/biocontainers/star:2.5.3a--0@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
  }

  output {
    String result = "STAR WORKES!"
  }

}
