workflow annotate_de_novo {
  File transcripts
  File blast_db
  File pfam
  Int threads
  String results_folder

  call transdecoder {
      input:
          transcripts = transcripts
  }

  call blastp as blast_orfs {
      input:
          query = transdecoder.orfs,
          db = blast_db,
          threads = threads
  }

  call identify_protein_domains {
    input:
        peptides = transdecoder.peptides,
        threads = threads,
        pfam = pfam
  }

  call copy as copy_transdecoder {
      input:
          files = [
            transdecoder.dir,
            transdecoder.peptides,
            transdecoder.cds,
            transdecoder.gff,
            transdecoder.positions,
            transdecoder.orfs
            ],
          destination = results_folder + "/transdecoder/"
  }

  call copy as copy_blast_orfs {
      input:
          files = [blast_orfs.out],
          destination = results_folder + "/blast/orfs/"
  }

  call copy as copy_domains {
      input:
          files = [identify_protein_domains.out],
          destination = results_folder + "/domains/"
  }

}

task transdecoder {

  File transcripts

  command {
    /opt/TransDecoder/TransDecoder.LongOrfs -t ${transcripts}
    /opt/TransDecoder/TransDecoder.Predict -t ${transcripts}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/transdecoder@sha256:5d2c702e7d430d8ea1c1dbb8d5706a0d4c208ace2c223f9fc43fb52fd0471f7c"
  }

  output {
    String name = basename(transcripts)
    File dir = name + ".transdecoder_dir"
    File peptides = name + ".transdecoder.pep"
    File cds = name + ".transdecoder.cds"
    File gff = name + ".transdecoder.gff3"
    File positions = name + ".transdecoder.bed"
    File orfs = dir + "/longest_orfs.pep"
  }
}

task blastp {

    File query
    File db
    Int threads

    command {
        makeblastdb -in ${db} -dbtype prot
        blastp -query ${query}  \
            -db ${db} \
            -num_threads ${threads} \
            -outfmt 6 > blastp.outfmt6
    }

    runtime {
        docker: "biocontainers/blast@sha256:71e8a6a9c7c697fd2fc9716e11a49b3a54267faafaaa4ad9ae32e5304e32834b"
    }

    output {
        File out = "blastp.outfmt6"
    }

}

task identify_protein_domains {

    File peptides
    File pfam
    Int threads

    command {
        hmmpress ${pfam}
        hmmscan --cpu ${threads} --domtblout hits.out ${pfam} ${peptides}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/hmmer@sha256:5b90f7d5c98d3adc41282159f69f013f7c02b6900ddd885ac1a654b19d704280"
    }

    output {
        File out = "hits.out"
    }

}

task signal_cleavage_prediction {

    File peptides

    command {
        signalp -f short -n signalp.out ${peptides}
    }

    runtime {
        #docker: "comics/hmmer@sha256:80d42008a2fb552f087c5320543f12610294dc174a28bf89b8141403dce19a33"
        docker: "quay.io/biocontainers/trinotate@sha256:070b1807fecb77c5997f6fe9c92e46a61824dd80e4f44a65958f00882f7c7e23"

        #TODO: find propert docker container
    }


    output {
        File out = "signalp.out"
    }

}

task transmembrate_regions_prediction {

    File proteins

    command {
        tmhmm --short ${proteins} tmhmm.out
    }

    runtime {
        docker: "quay.io/biocontainers/trinotate@sha256:070b1807fecb77c5997f6fe9c92e46a61824dd80e4f44a65958f00882f7c7e23"

        #TODO: find propert docker container
    }

    output {
        File out = "tmhmm.out"
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
