workflow downloads {

  String series
  Map[String, Array[String]] conditions

  scatter (condition in conditions) {
    call download{ input: series = series, condition = condition.left, samples = condition.right}
  }

}

task download {

  String series
  String condition
  Array[String] samples

  command {
    /opt/FastQC/fastqc ${file} -o .
  }

  output {
    #File out = sub(file, "\\.fastq.gz", "_fastqc.gz")
    File out = "${sampleName}_fastqc.zip"
  }
}
 #folder to download data to
