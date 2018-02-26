workflow vsearch {

  Int threads = 0
  File db
  File query
  String result_name
  String results_folder
  Float identity = 0.65

  call global_search {
      input:
        threads = threads,
        database = db,
        name = result_name,
        query = query,
        identity = identity
    }

  call copy as copy_results {
    input:
        files = [global_search.out],
        destination = results_folder
  }

  output {
       File out = copy_results.out[0]
  }

}

task global_search {

    File query
    File database
    String name
    Int threads
    Float identity

    command {
     vsearch --usearch_global ${query} --db ${database} --blast6out ${name}.blast6  --threads ${threads} --id ${identity}
    }
     #vsearch --usearch_global ${query} --db ${database} --threads ${threads} --id ${identity} --alnout alnout.txt


  runtime {
    docker: "quay.io/comp-bio-aging/vsearch:latest"
  }

  output {
       File out = name + ".blast6"
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