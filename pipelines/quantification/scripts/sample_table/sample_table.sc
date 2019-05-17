#! /usr/local/bin/amm

import $ivy.`com.github.pathikrit::better-files:3.4.0`
import $ivy.`com.nrinaudo::kantan.csv:0.5.0`
import $ivy.`com.nrinaudo::kantan.csv-generic:0.5.0`
import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._

implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

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
      print(s"GSE = " +series.name)
      //(series.name, series.children.filter(_.isDirectory))
      val samples = series.children.filter(s => s.isDirectory && s.nonEmpty).toList
      samples.flatMap {
        case gsm if gsm.isDirectory =>
          print("GSM = " + gsm.name)
          gsm.children.collect {
          case run if run.isDirectory && run.children.exists(_.name.contains("transcripts_abundance.tsv")) =>
            print(s"RUN = " +run.name)
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