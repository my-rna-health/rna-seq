workflow busco_check {

    String sequence_file
    File lineage
    Int threads
    String mode
    String name


    call busco {
        input:
           sequence = sequence_file,
           lineage = lineage,
           mode = mode,
           threads = threads,
           name = name
    }
}

task busco {

  File sequence
  String lineage
  String name
  String mode #--mode sets the assessment MODE: genome, proteins, transcriptome
  Int threads

  #python BUSCO.py -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]

  command {
    run_BUSCO.py -i ${sequence} -l /opt/busco/${lineage} -o ${name} -m ${mode} -c ${threads} -f
  }

  runtime {
    docker: "quay.io/comp-bio-aging/busco@sha256:fad688e8ddc8de05d75314d0278bc6fd51910466b9e96a6a494ff99fff1ac75b"
  }

  output {
    #File summary = "short_summary_" + name + ".txt"
    #File table = "full_table_" + name + ".tsv"
    #File missing_list = "missing_buscos_list_" + name + ".tsv"
    #File hmmer = "hmmer_output"
    #File blast = "blast_output"
    #File augustus = "augustus_output"
    #File single_copy_sequences = "single_copy_busco_sequences"
  }
}