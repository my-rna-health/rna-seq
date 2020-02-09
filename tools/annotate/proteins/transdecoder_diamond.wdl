workflow transdecoder_diamond {

  File transcripts
  File diamond_db
  Int threads
  String orfs_name
  String results_folder


  call transdecoder_orfs {
      input:
          transcripts = transcripts
  }

  call diamond_blast {
      input:
          query = transdecoder_orfs.orfs,
          database = diamond_db,
          threads = threads,
          name = orfs_name
  }

  call copy as copy_diamond_orfs {
      input:
          files = [diamond_blast.out],
          destination = results_folder + "/diamond/orfs/" + orfs_name + "/"
  }

  call transdecoder_predict {
      input:
        transcripts = transcripts,
        diamond_hits = diamond_blast.out
  }

  call copy as copy_transdecoder {
      input:
          files = [
            transdecoder_predict.dir,
            transdecoder_predict.peptides,
            transdecoder_predict.cds,
            transdecoder_predict.gff,
            transdecoder_predict.positions,
            transdecoder_predict.orfs
            ],
          destination = results_folder + "/transdecoder/" + orfs_name + "/"
  }

}

task transdecoder_orfs {

  File transcripts

  command {
    /opt/TransDecoder/TransDecoder.LongOrfs -t ${transcripts}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/transdecoder@sha256:5d2c702e7d430d8ea1c1dbb8d5706a0d4c208ace2c223f9fc43fb52fd0471f7c"
  }

  output {
    String name = basename(transcripts)
    File dir = name + ".transdecoder_dir"
    File orfs = dir + "/longest_orfs.pep"
    File cds = dir + "/longest_orfs.cds"
    File gff = dir + "/longest_orfs.gff3"
  }
}

task diamond_blast {

  Int threads
  File database
  File query
  String name

    command {
        diamond blastp -d ${database}  -q ${query} \
          --more-sensitive -o ${name}.m8 \
          -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
     }

  runtime {
    docker: "quay.io/comp-bio-aging/diamond:latest"
  }

  output {
       File out = name + ".m8"
  }

}


task transdecoder_predict {

  File transcripts
  File diamond_hits

  command {
    /opt/TransDecoder/TransDecoder.LongOrfs -t ${transcripts}
    /opt/TransDecoder/TransDecoder.Predict -t ${transcripts} --retain_blastp_hits ${diamond_hits}
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
