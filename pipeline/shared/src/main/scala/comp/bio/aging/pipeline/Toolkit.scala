package comp.bio.aging.pipeline


import cats.MonadCombine
import fr.hmil.roshttp.HttpRequest
import fr.hmil.roshttp.body.JSONBody.JSONObject

import scala.concurrent.ExecutionContext
import monix.execution.Scheduler.Implicits.global

import scala.util.{Failure, Success, Try}
import fr.hmil.roshttp.response.SimpleHttpResponse

import scala.concurrent.{Await, Future}
import scala.concurrent.duration._
import fr.hmil.roshttp.body.Implicits._
import fr.hmil.roshttp.body.JSONBody._
import fr.hmil.roshttp.body.{BodyPart, JSONBody, MultiPartBody, PlainTextBody}
import io.circe.Decoder.Result
import io.circe.generic.JsonCodec
import io.circe._
import io.circe.parser._
import io.circe.syntax._
import io.circe._
import io.circe.generic.semiauto._
import cats.implicits._

class ToolkitSRA(val base: String = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/") {


  def accession(id: String) = {
    val nums: String = id.replace("GSE", "")
    val ns: String = nums.take(nums.length - 3)
    base + s"GSE${ns}${"n" * 3}/"+id
    //GSE69263
  }

  /**
    *
    matrix/
    miniml/
    soft/
    suppl/
    */

}
