workflow quality_de_novo_genome {

  Int threads
  Int min_len
  Int q
  File reads_1
  File reads_2
  String results_folder #will be created if needed
  String adapter = "TruSeq3-PE"

  call report as initial_report_1 {
      input:
          sampleName = basename(reads_1, ".fastq.gz"),
          file = reads_1
  }

  call report as initial_report_2 {
      input:
          sampleName = basename(reads_2, ".fastq.gz"),
          file = reads_2
  }


  call trimmomatics {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        min_len = min_len,
        q = q,
        threads = threads,
        adapter = adapter
  }

  call report as report_trimmomatics_1 {
      input:
        sampleName = basename(trimmomatics.out1, ".fastq.gz"),
        file = trimmomatics.out1
      }

  call report as report_trimmomatics_2 {
      input:
        sampleName = basename(trimmomatics.out2, ".fastq.gz"),
        file = trimmomatics.out2
      }

  call copy as copy_trimmed {
    input:
        files = [trimmomatics.out1, trimmomatics.out2],
        destination = results_folder + "/trimmed/"
  }

  call copy as copy_initial_quality_reports {
    input:
        files = [initial_report_1.out, initial_report_2.out],
        destination = results_folder + "/quality/initial/"
  }

  call copy as copy_cleaned_quality_reports {
    input:
        files = [report_trimmomatics_1.out, report_trimmomatics_2.out],
        destination = results_folder + "/quality/cleaned/"
  }

  call multi_report {
    input:
        last_reports = copy_cleaned_quality_reports.out,
        folder = results_folder,
        report = "reports"
  }

  call copy as copy_multi_report {
      input:
          files = [multi_report.out],
          destination = results_folder
    }


  output {
    Array[File] out = copy_multi_report.out
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
  }

  output {
    File out = sampleName+"_fastqc.zip"
  }
}

task multi_report {

   File folder
   String report
   Array[File] last_reports #just a hack to make it wait for the folder to be created

   command {
        multiqc ${folder} --outdir ${report}
   }

   runtime {
        docker: "quay.io/comp-bio-aging/multiqc@sha256:20a0ff6dabf2f9174b84c4a26878fff5b060896a914d433be5c14a10ecf54ba3"
   }

   output {
        File out = report
   }
}


task trimmomatics {

    File reads_1
    File reads_2
    Int q
    Int min_len
    Int threads
    String adapter


    command {
       /usr/local/bin/trimmomatic PE \
            ${reads_1} \
            ${reads_2} \
            ${basename(reads_1, ".fastq.gz")}_trimmed.fastq.gz \
            ${basename(reads_1, ".fastq.gz")}_trimmed_unpaired.fastq.gz \
            ${basename(reads_2, ".fastq.gz")}_trimmed.fastq.gz \
            ${basename(reads_2, ".fastq.gz")}_trimmed_unpaired.fastq.gz \
            -threads ${threads} \
            ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/${adapter}.fa:2:30:10:1:TRUE \
            SLIDINGWINDOW:4:${q} MINLEN:${min_len}
    }

    runtime {
        docker: "quay.io/biocontainers/trimmomatic@sha256:bf4b0b2d2747670deeb9a6adc4a50aa923b830f0b02be47e82d1b848e1368078"
    }

    output {
        File out1 = basename(reads_1, ".fastq.gz") + "_trimmed.fastq.gz"
        File out1_unpaired = basename(reads_1, ".fastq.gz") + "_trimmed_unpaired.fastq.gz"
        File out2 = basename(reads_2, ".fastq.gz") + "_trimmed.fastq.gz"
        File out2_unpaired = basename(reads_2, ".fastq.gz") + "_trimmed_unpaired.fastq.gz"
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