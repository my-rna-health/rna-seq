package ro.biochim.aging.runner

import java.io.{File => JFile}

import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._

object StarAligner2 extends scala.App{

  //val host = "10.8.0.1" //localhost
  val host  = "pipelines.westeurope.cloudapp.azure.com"//"10.8.0.1"//"localhost"
  val port = "8000" //"38000"
  val urlEngine = s"http://${host}:${port}"
  val urlWorfklows = s"http://${host}:${port}/api"

  val engine = new CromwellClient(urlEngine, "v1")
  val client = new CromwellClient(urlWorfklows, "v1")

  val stats = client.waitFor(engine.getStats)


  val base = "/home/antonkulaga/rna-seq/workflows"
  val sourcePath = s"${base}/de-novo/assembly"
  //val workflow = s"${sourcePath}/star_aligner.wdl"
  val workflow = s"${sourcePath}/star_aligner2.wdl"
  //val inputs = s"${sourcePath}/inputs/2/kidney_star_ours.json"
  val inputs = s"${sourcePath}/inputs/2/liver_star_ours.json"
  //val inputs = s"${sourcePath}/inputs/their_index.json"

  val file = File(workflow)
  val input = File(inputs)

  def runWorkflow(): Status = {
    client.waitFor(client.postWorkflowFiles(file, input))
  }

  val status = runWorkflow()
  pprint.pprintln(status)
  println("OUTPUTS")
  val outputs = client.waitFor(client.getAllOutputs())
  pprint.pprintln(outputs)

}