#! /usr/local/bin/amm
import coursier.maven.MavenRepository
interp.repositories() ++= Seq(
  MavenRepository("https://dl.bintray.com/comp-bio-aging/main/")
)
@
import $ivy.`com.github.pathikrit::better-files:3.4.0`
import $ivy.`com.nrinaudo::kantan.csv:0.5.0`
import $ivy.`com.nrinaudo::kantan.csv-generic:0.5.0`
import $ivy.`group.aging-research::geo-fetch:0.0.3`
import $ivy.`io.circe::circe-optics:0.11.0`
import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._

implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

lazy val key ="0a1d74f32382b8a154acacc3a024bdce3709"

case class Run(gse: String, gsm: String, run: String, genes: String, transcripts: String, quant: String, organism: String)


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
  val runs = rt.children.filter(_.isDirectory).flatMap {
    series =>
      println("============")
      println(s"GSE = " +series.name)
      //(series.name, series.children.filter(_.isDirectory))
      val samples = series.children.filter(s => s.isDirectory && s.nonEmpty).toList
      samples.flatMap {
        case gsm if gsm.isDirectory =>
          println("GSM = " + gsm.name)
          if(! (gsm.children.exists(_.name == gsm.name + ".json") && gsm.children.exists(_.name == gsm.name + "runs.tsv"))  )
          {
            println("cannot find " + gsm.name + ".json, creating it from scratch!")
            geo.cli.MainCommand.fetchGSM(gsm.name, key, gsm.pathAsString + "/" + gsm.name + ".json", gsm.pathAsString + "/" + gsm.name +"_runs.tsv", true)
          }

          gsm.children.collect {
            case run if run.isDirectory && run.children.exists(_.name.contains("transcripts_abundance.tsv")) =>
              println(series.name + "/" + gsm.name + "/" + run.name)
              Run(
                series.name,
                gsm.name,
                run.name,
                getPathIf(run)(_.name.contains("genes_abundance.tsv")),
                getPathIf(run)(_.name.contains("transcripts_abundance.tsv")),
                getPathIf(run)(_.name.contains("quant_")),
                ""
              )
          }
      }
  }.toList
  index.toIO.asCsvWriter[Run](config.withHeader(true)).write(runs)
}