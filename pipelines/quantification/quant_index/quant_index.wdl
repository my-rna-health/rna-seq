version development

struct Gentrome {
    String species
    File genome
    File transcriptome
    String version
    String? subversion
}

workflow quant_index {
    input {
        Array[Gentrome] references
        String indexes_folder
        Int threads = 9
        String max_memory = "24G"
    }

    scatter (gentrome in references) {
        String organism = sub(gentrome.species, " ", "_")

        call make_decoys {
           input: reference = gentrome
        }

        String full_index_name = gentrome.version + if(defined(gentrome.subversion)) then "_" + gentrome.subversion else ""

        call salmon_index  {
           input:
               indexName =  full_index_name,
               gentrome = make_decoys.gentrome,
               decoys = make_decoys.decoys,
               p = threads,
               max_memory = max_memory
        }

        call copy_folder {
         input:
             dir = salmon_index.out,
             destination = indexes_folder + "/" + organism
        }
    }

    output {
        Array[File] indexes = copy_folder.out
    }

}

task make_decoys {
    input {
       Gentrome reference
    }

    String name = basename(reference.genome)

    command {
        grep "^>" ~{reference.genome} | cut -d " " -f 1 > ~{name}_decoys.txt
        sed -i -e 's/>//g' ~{name}_decoys.txt
        cat ~{reference.transcriptome} ~{reference.genome} > ~{name}_gentrome.fa
    }

    output {
        File decoys = name + "_decoys.txt"
        File gentrome = name + "_gentrome.fa"
    }
}

task salmon_index {
    input {
        File gentrome
        File decoys
        String indexName
        Int p = 3
        String max_memory= "48G"
    }

  command {
    salmon index -t ~{gentrome} -d ~{decoys} -p ~{p} -i ~{indexName}
  }

  runtime {
    docker: "quay.io/biocontainers/salmon@sha256:0aea466ba3eae62cb3ea8077e30b6212ca152f73e9520ca19f7421c5be519ef9" #1.1.0--hf69c8f4_0
    maxRetries: 3
    docker_memory: "${max_memory}"
    docker_cpu: "${p}"
  }

  output {
    Directory out = indexName
  }

}

task copy_folder {
    input {
        Directory dir
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{dir} ~{where}
    }

    output {
        File out = where
    }
}
