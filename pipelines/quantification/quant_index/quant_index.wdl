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
               decoys = make_decoys.decoys
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
        grep "^>" <(zcat ~{reference.genome}) | cut -d " " -f 1 > ~{name}_decoys.txt
        sed -i -e 's/>//g' ~{name}_decoys.txt
        cat ~{reference.transcriptome} ~{reference.genome} > ~{name}_gentrome.fa.gz
    }

    output {
        File decoys = name + "_decoys.txt"
        File gentrome = name + "_gentrome.fa.gz"
    }
}

task salmon_index {
    input {
        File gentrome
        File decoys
        String indexName
        Int p = 12
    }

  command {
    salmon index -t ~{gentrome} -d ~{decoys} -p ~{p} -i ~{indexName}
  }

  runtime {
    docker: "quay.io/biocontainers/salmon@sha256:7182223f62fad3c1049342cc686f1c5b6991c6feaa7044ed3532dbbd4e126533" #1.0.0--hf69c8f4_0
    maxRetries: 3
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