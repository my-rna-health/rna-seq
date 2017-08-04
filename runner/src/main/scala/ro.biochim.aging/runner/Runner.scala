package ro.biochim.runner


import java.io.{File => JFile}
import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._

object Runner extends scala.App{

  val port = "38000"
  //val client = CromwellClient.localhost
  val client = new CromwellClient(s"http://localhost:${port}/api", "v1")

  import fr.hmil.roshttp.body.Implicits._


  val stats = client.waitFor(client.getStats)
  val base = "/home/antonkulaga/rna-seq/workflows"
  val sourcePath = s"${base}/rna-seq"
  val workflow = s"${sourcePath}/main_workflow.wdl"
  val inputs = s"${sourcePath}/inputs/worm.input.json"
  val subs = s"${sourcePath}/subs"

  val file = File(workflow)

  val input = File(inputs)

  def runWorkflow(): Status = {
    client.waitFor(client.postWorkflowFiles(file, input, File(subs)))
  }

  val status = runWorkflow()
  pprint.pprintln(status)
  println("OUTPUTS")
  val outputs = client.waitFor(client.getAllOutputs())
  pprint.pprintln(outputs)
}
