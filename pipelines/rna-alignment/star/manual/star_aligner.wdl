version development

workflow StarAligner {

  input {
    File index_dir
    File reads_1
    File reads_2
    String results_folder #will be created if needed
    Float threshold  = 0.66
    Int threads = 8
  }


  call star_align {
      input:
        reads_1 = reads_1,
        reads_2 = reads_2,
        index_dir = index_dir,
        threads = threads,
        threshold = threshold
    }

  call copy as copy_results {
    input:
        files = [star_align.log, star_align.out, star_align.junctions],
        destination = results_folder
  }

  output {
    Array[File] results = copy_results.out
  }

}

task star_align {
    input {
        File reads_1
        File reads_2
        File index_dir
        Int  threads
        Float threshold
    }


  command {
    /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --genomeDir ~{index_dir} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand gunzip -c \
        --outFilterMatchNminOverLread ~{threshold} \
        --readFilesIn ~{reads_1} ~{reads_2}
  } # --outSAMtype BAM SortedByCoordinate

  runtime {
    docker: "quay.io/biocontainers/star@sha256:f9b0406354ff2e5ccfadaef6fde6367c7bcb4bdc7e67920f0f827a6ff6bf4fb5" #2.7.2b--0
  }

  output {
    File out = "Aligned.sortedByCoord.out.bam" #"Aligned.out.sam"
    File log = "Log.final.out"
    File junctions = "SJ.out.tab"
  }

}

task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{where}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}