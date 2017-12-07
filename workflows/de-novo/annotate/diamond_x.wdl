workflow Diamond_X {

  Int threads
  File db
  String results_name
  String results_folder
  String db_name = basename(db, ".m8")





  call diamond_blastx {
      input:
        threads = threads,
        peptides = db,
        db_name = db_name
    }

  call copy as copy_results {
    input:
        files = [diamond_blastp.out],
        destination = results_folder
  }

}

task diamond_blastx {

  Int threads
  File orfs
  File db_name
  String results_name

    command {
        diamond blastx -d ${db_name} -q ${orfs} --more-sensitive -o ${results_name}.m8
     }

  runtime {
    docker: "quay.io/biocontainers/diamond@sha256:ad7ed429a1a0ee95e88c29b10b44032ce1ab23c9ef91bf49e9062aa10ec91231"
  }

  output {
       File out = results_name + ".m8"
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