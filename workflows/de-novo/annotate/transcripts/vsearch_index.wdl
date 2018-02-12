workflow vsearch_index {

  File db
  String name
  String indexes_folder = "/pipelines/indexes/vsearch"

  call vsearch_make_index {
    input:
      fasta = db,
      name = name
  }

  call copy as copy_results {
    input:
        files = [vsearch_make_index.out],
        destination = indexes_folder
  }

  output {
    Array[File] results = copy_results.out
  }

}

task vsearch_make_index {

  File fasta
  String name

    command {
        vsearch --makeudb_usearch ${fasta} --output ${name}.udb
     }

  runtime {
    docker: "quay.io/comp-bio-aging/vsearch:latest"
  }

  output {
       File out = name + ".udb"
       String db_name = name
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