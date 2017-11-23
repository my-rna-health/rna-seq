workflow Diamond {

  Int threads
  File db
  String name
  String results_folder

  call diamond_index {
      input:
        threads = threads,
        fasta = db,
        name = name
    }

  call diamond_blastp {
      input:
        threads = threads,
        peptides = db,
        name = name
    }

  call copy as copy_results {
    input:
        files = [diamond_blastp.out],
        destination = results_folder
  }

  output {
    Array[File] results = copy_results.out
  }

}

task diamond_index {

  Int threads
  File fasta
  String name

    command {
        diamond \
        --threads ${threads}\
        makedb  --in ${fasta} ${name}
     }

  runtime {
    docker: "quay.io/biocontainers/diamond@sha256:ad7ed429a1a0ee95e88c29b10b44032ce1ab23c9ef91bf49e9062aa10ec91231"
  }

  output {
       File out = name + ".m8"
       String db_name = name
  }

}

task diamond_blastp {

  Int threads
  File peptides
  File name

    command {
        diamond blastp -d ${name} -q ${peptides} -o matches.m8
     }

  runtime {
    docker: "quay.io/biocontainers/diamond@sha256:ad7ed429a1a0ee95e88c29b10b44032ce1ab23c9ef91bf49e9062aa10ec91231"
  }

  output {
       File out = "matches.m8"
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