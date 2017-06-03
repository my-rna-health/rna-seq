task star {

  Int  numberOfThreads = 8
  String sampleName
  File file
  File genomeDir

  command {
    /usr/local/bin/STAR --runThreadN ${numberOfThreads} --genomeDir ${genomeDir} --readFilesCommand gunzip -c --readFilesIn ${file} --quantMode TranscriptomeSAM
  }

  runtime {
    docker: "quay.io/biocontainers/star@sha256:352f627075e436016ea2c38733b5c0096bb841e2fadcbbd3d4ae8daf03ccdf1b"
  }

  output {
    Map[String, File] out = {"sample": file, "file": "Aligned.out.sam", "log": "Log.final.out"}
  }

}
