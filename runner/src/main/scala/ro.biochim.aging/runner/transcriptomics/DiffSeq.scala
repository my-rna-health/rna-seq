package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import ro.biochim.aging.runner.core.BasicRunner

object DiffSeq extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "rna-seq","main_workflow.wdl"){
  runWithSubs("worm.input.json")
}
