package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import ro.biochim.aging.runner.core.BasicRunner

object StarAligner extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "de-novo/star","star_aligner.wdl"){
  //run("liver_star_ours.json")
  //run("liver_star_theirs.json")
  //run("kidney_star_theirs.json")
  //run("liver_star_ours.json")
  //run("kidney_star_bowhead_whale.json")
}
