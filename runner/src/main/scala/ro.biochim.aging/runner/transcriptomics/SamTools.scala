package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import ro.biochim.aging.runner.core.BasicRunner

object SamTools extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "de-novo/samtools","samtools.wdl"){
  run("merge_ours.json")

}
