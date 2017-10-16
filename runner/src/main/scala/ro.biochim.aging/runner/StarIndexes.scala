package ro.biochim.aging.runner

import java.io.{File => JFile}

import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._

object StarIndexes extends scala.App{

  //val host = "10.8.0.1" //localhost
  val host  = "13.81.207.19"//"pipelines.westeurope.cloudapp.azure.com"//"10.8.0.1"//"localhost"
  val port = "8000" //"38000"
  val urlEngine = s"http://${host}:${port}"
  val urlWorfklows = s"http://${host}:${port}/api"

  val engine = new CromwellClient(urlEngine, "v1")
  val client = new CromwellClient(urlWorfklows, "v1")

  import fr.hmil.roshttp.body.Implicits._


  //println("server is: "+url)
  //val stats = client.waitFor(engine.getStats)

  val base = "/home/antonkulaga/rna-seq/workflows"
  val sourcePath = s"${base}/indexes/star"
  //val workflow = s"${sourcePath}/star_aligner.wdl"
  val workflow = s"${sourcePath}/star_index.wdl"
  //val inputs = s"${sourcePath}/inputs/our_index.json"
  //val inputs = s"${sourcePath}/inputs/their_index.json"
  //val inputs = s"${sourcePath}/inputs/bowhead_index.json"
  val inputs = s"${sourcePath}/inputs/human_index.json"


  //val subs = s"${sourcePath}/subs"
  val stats = client.waitFor(engine.getStats)
  println(stats)

  val file = File(workflow)
  val input = File(inputs)

  def runWorkflow(): Status = {
    val wdl = file.lines.mkString("\n")
    val i= if(input.exists) input.lines.mkString("\n") else ""
    println(wdl)
    println(i)
    client.waitFor(client.postWorkflowFiles(file, input))
  }



  val status = runWorkflow()
  pprint.pprintln(status)
  println("OUTPUTS")
  val outputs = client.waitFor(client.getAllOutputs())
  pprint.pprintln(outputs)


}
