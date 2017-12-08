workflow Diamond_Blast {

  Int threads
  File db
  String result_name
  String results_folder
  String mode # = "blastx"

  call diamond_blastp {
      input:
        threads = threads,
        database = db,
        name = result_name
    }

  call copy as copy_results {
    input:
        files = [diamond_blastp.out],
        destination = results_folder
  }

  output {
       File out = copy_results[0]
  }

}

task diamond_blast {

  Int threads
  File database
  String name
  String mode

    command {
        diamond ${mode} -d ${database} -q ${database} \
          --more-sensitive -o ${name}.m8 \
          -f 6 id sseqid qseq score pident staxids stitle qcovhsp qtitle \
     }

  runtime {
    docker: "quay.io/comp-bio-aging/diamond:latest"
  }

  output {
       File out = name + ".m8"
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