package ro.biochim.aging.runner.genomics

import java.io.{File => JFile}

import ro.biochim.aging.runner.core.BasicRunner
import ro.biochim.aging.runner.transcriptomics.SamTools.run

object ABYSS extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "de-novo/genome-assembly","abyss.wdl"){
  run("ours.json")
}