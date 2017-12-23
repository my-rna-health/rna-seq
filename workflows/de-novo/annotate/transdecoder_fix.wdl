workflow transdecoder_fix {

  File transcripts
  File pfam_hits
  File diamond_hits

  String results_folder

  call transdecoder {
      input:
          transcripts = transcripts,
          pfam_hits = pfam_hits,
          diamond_hits = diamond_hits

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

}

task transdecoder {

  File transcripts
  File pfam_hits
  File diamond_hits

  command {
    /opt/TransDecoder/TransDecoder.LongOrfs -t ${transcripts}
    /opt/TransDecoder/TransDecoder.Predict -t ${transcripts} --retain_pfam_hits ${pfam_hits} --retain_blastp_hits ${diamond_hits}
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
