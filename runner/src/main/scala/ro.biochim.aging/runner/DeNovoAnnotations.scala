package ro.biochim.aging.runner

import java.io.{File => JFile}
import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._

object DeNovoAnnotations extends scala.App{

  //val host = "10.8.0.1" //localhost
  val host  = "10.8.0.1"//"localhost"
  val port = "8000" //"38000"
  val url = s"http://${host}:${port}/api"

  val engine = new CromwellClient(s"http://${host}:${port}", "v1")
  val client = new CromwellClient(url, "v1")

  import fr.hmil.roshttp.body.Implicits._


  //println("server is: "+url)
  val stats = client.waitFor(engine.getStats)

  val base = "/home/antonkulaga/rna-seq/workflows"
  val sourcePath = s"${base}/de-novo"
  val workflow = s"${sourcePath}/annotations_de_novo.wdl"
  val inputs = s"${sourcePath}/inputs/annotations_de_novo.json"
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
