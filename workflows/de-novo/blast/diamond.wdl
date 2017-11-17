workflow Diamond {

  File index_dir
  Int threads

  call diamond_index {
      input:
        threads = threads
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

task diamond_index {

  Int threads

  command {
        diamond \
        --threads ${threads}\
        makedb
     } # --outSAMtype BAM SortedByCoordinate

  runtime {
    docker: "quay.io/biocontainers/diamond@sha256:ad7ed429a1a0ee95e88c29b10b44032ce1ab23c9ef91bf49e9062aa10ec91231"
  }

  output {

  }

}