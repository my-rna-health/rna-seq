workflow annotate_de_novo {
  File transcripts
  File blast_db
  Int threads

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

}

task transdecoder {

  File transcripts

  command {
    /opt/TransDecoder/TransDecoder.LongOrfs -t ${transcripts}
    /opt/TransDecoder/TransDecoder.Predict -t ${transcripts}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/transdecoder@sha256:82e9f9372c0624a197d0fc97f477412ee46e0b20edb211fa604d074975eec46a"
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

task trinnotate_index {


 runtime {
    docker: "quay.io/biocontainers/trinotate@sha256:070b1807fecb77c5997f6fe9c92e46a61824dd80e4f44a65958f00882f7c7e23"
 }

 output {
    File db = "Trinotate.sqlite"
    File uniprot = "uniprot_sprot.pep"
    File Pfam = "Pfam-A.hmm.gz"
 }

}

task identify_protein_domains {

    File predictions
    Int threads

    command {
        hmmscan --cpu ${threads} --domtblout TrinotatePFAM.out Pfam-A.hmm transdecoder.pep > pfam.log
    }

    runtime {
        #docker: "comics/hmmer@sha256:80d42008a2fb552f087c5320543f12610294dc174a28bf89b8141403dce19a33"
        docker: "quay.io/biocontainers/trinotate@sha256:070b1807fecb77c5997f6fe9c92e46a61824dd80e4f44a65958f00882f7c7e23"
    }

    output {

    }

}

task signal_cleavage_prediction {

    command {
        signalp -f short -n signalp.out transdecoder.pep
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
