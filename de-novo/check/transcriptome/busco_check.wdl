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
    docker: "quay.io/comp-bio-aging/busco@sha256:464c127a0eb3ad074a72ae0f4841d86e295de15b9cff71dc1fd722976ddab9da"
  }

  output {

    String prefix = "run_"+name + "/"
    File summary = prefix + "short_summary_" + name + ".txt"
    File table = prefix + "full_table_" + name + ".tsv"
    File missing_list = prefix + "missing_busco_list_" + name + ".tsv"
    File hmmer = prefix + "hmmer_output"
    File blast = prefix + "blast_output"
    File proteins = prefix + "translated_proteins"
    #File augustus = "augustus_output"
    #File single_copy_sequences = "single_copy_busco_sequences"
  }
}