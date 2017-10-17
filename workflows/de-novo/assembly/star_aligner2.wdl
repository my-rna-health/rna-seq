#TODO: merge with StarAligner
workflow StarAligner2 {

  String title
  File index_dir
  Int threads
  File reads_1
  File reads_2
  Array[File] junctions

  call star_align2 {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        index_dir = index_dir,
        threads = threads,
        junctions_1 = junctions
    }

  call copy as copy_results {
    input:
        files = [star_align2.log, star_align2.out, star_align2.junctions],
        destination = results_folder
  }

  output {
    Array[File] results = copy_results.out
    File name = title
  }

}

task star_align2 {

  File reads_1
  File reads_2
  File index_dir
  Int  threads
  Array[File] junctions_1


  command {
    /usr/local/bin/STAR \
        --runThreadN ${threads} \
        --genomeDir ${index_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${reads_1} ${reads_2} \
        --sjdbFileChrStartEnd ${sep=' ' junctions_1}
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

task copy {
    Array[File] files
    String destination

    command {
        mkdir -p ${destination}
        cp -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}
