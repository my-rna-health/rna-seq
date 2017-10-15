workflow StarAligner {

  String title
  File index_dir
  Int threads
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
        min_len = 36,
        q = 19,
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

  call star_align {
      input:
        reads_1 = trimmomatics.out1,
        reads_2 = trimmomatics.out2,
        index_dir = index_dir,
        threads = threads
    }

  call copy as copy_results {
    input:
        files = [report_trimmomatics_1.out, report_trimmomatics_2.out, star_align.log, star_align.out, star_align.junctions],
        destination = results_folder
  }

  output {
    Array[File] results = copy_results.out
    File name = title
  }

}

task star_align {

  File reads_1
  File reads_2
  File index_dir
  Int  threads


  command {
    /usr/local/bin/STAR \
        --runThreadN ${threads} \
        --genomeDir ${index_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${reads_1} ${reads_2}
  } # --outSAMtype BAM SortedByCoordinate

  runtime {
    docker: "quay.io/biocontainers/star@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
  }

  output {
    File out = "Aligned.out.sam" #"Aligned.sortedByCoord.out.bam"
    File log = "Log.final.out"
    File junctions = "SJ.out.tab"
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
    File destination

    command {
        mkdir -p ${destination}
        cp -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}
