workflow downloads {

  File folder #folder to download data to
  Map[String, Array[String]] conditions

}

task report {

  String sampleName
  File file

  command {
    /opt/FastQC/fastqc ${file} -o .
  }

  runtime {
    docker: "quay.io/ucsc_cgl/fastqc@sha256:86d82e95a8e1bff48d95daf94ad1190d9c38283c8c5ad848b4a498f19ca94bfa"
    #docker: "quay.io/biocontainers/fastqc@sha256:bb57a4deeec90633e746afbc38c36fdb202599fe71f9557b94652e9c8f3c1a02"
  }

  output {
    #File out = sub(file, "\\.fastq.gz", "_fastqc.gz")
    File out = "${sampleName}_fastqc.zip"
  }
}
 #folder to download data to
