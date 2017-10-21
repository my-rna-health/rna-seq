package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import comp.bio.aging.cromwell.client._

import java.io.{File => JFile}

import better.files._
import comp.bio.aging.cromwell.client._
import ro.biochim.aging.runner.core.BasicRunner

object SalmonIndexes extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "de-novo/indexes/salmon","indexes.wdl"){

  //run("fly.json")
  //run("mouse.json")
  //run("human.json")
  run("worm.json")
}
