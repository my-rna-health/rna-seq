package ro.biochim.aging.runner.transcriptomics

import java.io.{File => JFile}

import better.files._
import comp.bio.aging.cromwell.client._
import ro.biochim.aging.runner.core.BasicRunner
import ro.biochim.aging.runner.transcriptomics.DeNovoAnnotations.run

object Trinity extends BasicRunner(
  "/home/antonkulaga/rna-seq/workflows",
  "de-novo/trinity","assembly_trinity.wdl"){

  run("genome_guided.json")


}
