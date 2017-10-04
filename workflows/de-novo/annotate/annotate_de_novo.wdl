workflow annotate_de_novo {
  File transcripts

  call transdecoder {
      input:
          transcripts = transcripts
  }
}

task transdecoder {


  File transcripts


  command {
    /opt/TransDecoder/TransDecoder.LongOrfs -t ${transcripts}
    /opt/TransDecoder/TransDecoder.Predict -t ${transcripts}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/transdecoder:latest"
  }

  output {
    String name = basename(transcripts)
    File dir = name + ".transdecoder_dir"
    File peptides = name + ".transdecoder.pep"
    File cds = name + ".transdecoder.cds"
    File gff = name + ".transdecoder.gff3"
    File positions = name + ".transdecoder.bed"
  }
}
