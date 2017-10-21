package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import better.files._
import comp.bio.aging.cromwell.client._

import java.io.{File => JFile}

import ro.biochim.aging.runner.core.BasicRunner

object CheckBusco extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "de-novo/check","check_de_novo.wdl"){
  run("fly_metazoan.json")
}
