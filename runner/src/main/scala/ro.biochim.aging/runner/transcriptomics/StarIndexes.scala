package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import better.files._
import comp.bio.aging.cromwell.client._
import ro.biochim.aging.runner.core.BasicRunner

object StarIndexes extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "indexes/star","star_index.wdl"){
  run("liver_star_ours.json")
  //run("our_index.json")
  //run("their_index.json")
  //run("bowhead_index.json")
  //run("human_index.json")
  //run("cow_index.json")
}
