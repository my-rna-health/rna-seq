package ro.biochim.runner


import java.io.{File => JFile}
import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._


object Indexes extends scala.App{
  import better.files._
  //val host = "10.8.0.1" //localhost
  val host  = "10.8.0.1"//"localhost"
  val port = "8000" //"38000"

  val client = new CromwellClient(s"http://${host}:${port}/api", "v1")

  val stats = client.waitFor(client.getStats)
  val base = "/home/antonkulaga/rna-seq/workflows"
  val sourcePath = s"${base}/indexes"
  val workflow = s"${sourcePath}/indexes.wdl"
  val subs = s"${sourcePath}/subs"


  //val inputs = s"${sourcePath}/inputs/fly.json"
  //val inputs = s"${sourcePath}/inputs/mouse.json"
  //val inputs = s"${sourcePath}/inputs/human.json"
  val inputs = s"${sourcePath}/inputs/worm.json"





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
