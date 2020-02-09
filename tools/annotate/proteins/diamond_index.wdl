workflow Diamond_Index {

  File db
  String name
  String results_folder

  call diamond_index {
    input:
      fasta = db,
      name = name
  }

  call copy as copy_results {
    input:
        files = [diamond_index.out],
        destination = results_folder
  }

  output {
    Array[File] results = copy_results.out
  }

}

task diamond_index {

  File fasta
  String name

    command {
        diamond makedb --in ${fasta} -d ${name}
     }

  runtime {
    docker: "quay.io/comp-bio-aging/diamond:latest"
  }

  output {
       File out = name + ".dmnd"
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