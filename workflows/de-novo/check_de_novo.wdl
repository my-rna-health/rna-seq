workflow check_de_novo {

    File transcriptomeIndex
}

task quality {


}


task busco {


  File sequence
  String lineage
  String name
  String mode #--mode sets the assessment MODE: genome, proteins, transcriptome

  #python BUSCO.py -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]

  command {
    run_busco -i ${sequence} -l ${lineage} -o ${name} -m ${mode}
  }

  runtime {
    docker: "quay.io/biocontainers/busco:3.0.2--py35_2"

  output {
    File out = "${sampleName}_trimmed.fastq.gz"
  }
}
