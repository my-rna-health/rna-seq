workflow Diamond_Blast {

  Int threads
  File db
  File query
  String result_name
  String results_folder
  String mode = "blastx"

  call diamond_blast {
      input:
        threads = threads,
        database = db,
        name = result_name,
        query = query,
        mode = mode

    }

  call copy as copy_results {
    input:
        files = [diamond_blast.out],
        destination = results_folder
  }

  output {
       File out = copy_results.out[0]
  }

}

task diamond_blast {

  Int threads
  File database
  File query
  String name
  String mode

    command {
        diamond ${mode} -d ${database}  -q ${query} \
          --more-sensitive -o ${name}.m8 \
          -f 6 sseqid qseq score pident stitle qcovhsp qtitle \
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