package ro.biochim.aging.runner

import java.io.{File => JFile}
import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._

object DeNovoQuality extends scala.App{

  //val host = "10.8.0.1" //localhost
  val host  = "pipelines1.westeurope.cloudapp.azure.com"//"10.8.0.1"//"localhost"
  val port = "8000" //"38000"
  val urlEngine = s"http://${host}:${port}"
  val urlWorfklows = s"http://${host}:${port}/api"

  val engine = new CromwellClient(urlEngine, "v1")
  val client = new CromwellClient(urlWorfklows, "v1")

  val stats = client.waitFor(engine.getStats)


  val base = "/home/antonkulaga/rna-seq/workflows"
  val sourcePath = s"${base}/de-novo/quality"
  val workflow = s"${sourcePath}/quality_de_novo.wdl"
  val inputs = s"${sourcePath}/inputs/wilver.json"
  //val inputs = s"${sourcePath}/inputs/mother_kidney.json"
  //val subs = s"${sourcePath}/subs"

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
