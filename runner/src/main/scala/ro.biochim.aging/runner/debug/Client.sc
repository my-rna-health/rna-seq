

import java.io.{File => JFile}

import scala.concurrent.duration._
import better.files._
import comp.bio.aging.cromwell.client._
import fr.hmil.roshttp.body.JSONBody._

//val host = "10.8.0.1" //localhost
val host  = "13.81.207.19"//"pipelines.westeurope.cloudapp.azure.com"//"10.8.0.1"//"localhost"
val port = "8000" //"38000"
val urlEngine = s"http://${host}:${port}"
val urlWorfklows = s"http://${host}:${port}/api"

val engine = new CromwellClient(urlEngine, "v1")
val client = new CromwellClient(urlWorfklows, "v1")

val succeeded = client.waitFor(client.getQuery(WorkflowStatus.Succeeded))

pprint.pprintln(succeeded)
