#! /usr/local/bin/amm
import $exec.dependencies

import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._
import io.circe.optics.JsonPath.root
import io.circe._, io.circe.parser._
import scala.util._

import io.circe.optics.JsonPath._
val title = root.title.string
val organism = root.organism.name.string
val taxid = root.organism.taxid.string
val sequencer = root.sequencer.string
val strategy = root.library.strategy.string
val selection = root.library.selection.string
val source = root.library.source.string
val extraction_source   = root.extraction.source.string
val molecule = root.extraction.molecule.string
val protocol = root.extraction.protocol.string

implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

lazy val key ="0a1d74f32382b8a154acacc3a024bdce3709"
object Run {
  implicit val sampleCodec: HeaderCodec[Run] = HeaderCodec.caseCodec(
    "gse",	"gsm",	"run",	"organism",	"taxid",
    "title",	"sequencer",	"strategy",	"selection",
    "source",	"extraction_source", "molecule", "genes", "transcripts", "quant", "protocol")(Run.apply)(Run.unapply)
}
case class Run(gse: String, gsm: String, run: String, organism: String,
               taxid: String, title: String,  sequencer: String,
               strategy: String, selection: String,
               source: String, extraction_source: String,
               molecule: String, genes: String, transcripts: String, quant: String, protocol: String
              )
object SampleInfo{
  def empty(gsm: String, gse: String): SampleInfo = {
    SampleInfo(gsm, gse, "", "", "", "", "", "", "", "", "", "")
  }
}
case class SampleInfo(gsm: String, gse: String, title: String, organism: String, taxid: String,
                      sequencer: String, strategy: String, selection: String, source: String, extraction_source: String, molecule: String, protocol: String)

def getPathIf(dir: File)(fun: File => Boolean) = dir.children.collectFirst{ case f if fun(f) => f.pathAsString }.getOrElse("")

/**
  * Builds an index of quantified GSE/GSM/SRR-s
  * @param index index file name for samples
  * @param root root folder (/data/samples by default)
  * @return
  */
@main
def main(index: Path = Path("/data/samples/index.tsv"), root: Path = Path("/data/samples")) = {
  val fl = index.toIO.toScala
  val rt = root.toIO.toScala
  val runs: List[Run] = rt.children.filter(_.isDirectory).flatMap {
    series =>
      println("============")
      println(s"GSE = " +series.name)
      //(series.name, series.children.filter(_.isDirectory))
      val samples = series.children.filter(s => s.isDirectory && s.nonEmpty).toList
      samples.flatMap {
        case gsm if gsm.isDirectory && gsm.children.exists(f=> f.isDirectory &&
          f.children.exists(child=>child.name.contains("_transcripts_abundance.tsv"))
        ) =>
          println("GSM = " + gsm.name)
          if(! (gsm.children.exists(_.name == gsm.name + ".json") && gsm.children.exists(_.name == gsm.name + "_runs.tsv"))  )
          {
            println("cannot find " + gsm.name + ".json, creating it from scratch!")
            Try {
              geo.cli.MainCommand.fetchGSM(gsm.name, key, gsm.pathAsString + "/" + gsm.name + ".json", gsm.pathAsString + "/" + gsm.name + "_runs.tsv", true)
            } match {
              case Failure(th) => println("could not create GSM because of: " + th.toString)
               case _ => println(gsm.name + ".json" + " successfully created!")
            }
          }
          val sample = if(gsm.children.exists(f=> f.name == gsm.name + ".json" && !f.name.contains("PRJ"))){
            parse((gsm / (gsm.name + ".json")).contentAsString) match{
              case Left(io.circe.ParsingFailure(message: String, underlying: Throwable)) =>
                println(s"could not parse ${gsm.pathAsString} because: $message and ${underlying}")
                SampleInfo.empty(gsm.name, series.name)
              case Right(json) =>
                SampleInfo(
                  gsm.name,
                  series.name,
                  title.getOption(json).getOrElse(""),
                  organism.getOption(json).getOrElse(""),
                  taxid.getOption(json).getOrElse(""),
                  sequencer.getOption(json).getOrElse(""),
                  strategy.getOption(json).getOrElse(""),
                  selection.getOption(json).getOrElse(""),
                  source.getOption(json).getOrElse(""),
                  extraction_source.getOption(json).getOrElse(""),
                  molecule.getOption(json).getOrElse(""),
                  protocol.getOption(json).getOrElse("")
                )

            }
          }
          else SampleInfo.empty(gsm.name, series.name)

          gsm.children.collect {
            case run if run.isDirectory && run.children.exists(_.name.contains("_transcripts_abundance.tsv")) =>
              println(series.name + "/" + gsm.name + "/" + run.name)
              Run(
                series.name,
                gsm.name,
                run.name,
                sample.organism, sample.taxid,
                sample.title, sample.sequencer, sample.strategy,
                sample.selection, sample.source, sample.extraction_source, sample.molecule,
                getPathIf(run)(_.name.contains("genes_abundance.tsv")),
                getPathIf(run)(_.name.contains("transcripts_abundance.tsv")),
                getPathIf(run)(_.name.contains("quant_")),
                sample.protocol
              )
          }
        case gsm =>
          println(s"GSM ${gsm.name} does not seem to have SRR-s inside!")
          Nil
      }
  }.toList
  index.toIO.asCsvWriter[Run](config.withHeader).write(runs)
  println("INDEX SUCCESSFULLY CREATED at " + index.toIO.toScala.pathAsString)
}