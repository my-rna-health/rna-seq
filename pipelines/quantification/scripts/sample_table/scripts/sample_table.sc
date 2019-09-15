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
import better.files.File
import geo.cli._
import geo.fetch._


implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

lazy val key ="0a1d74f32382b8a154acacc3a024bdce3709"

def getPathIf(dir: File)(fun: File => Boolean) = dir.children.collectFirst{ case f if fun(f) => f.pathAsString }.getOrElse("")



/**
  * Write indexes with runs to the disck
  * @param index
  * @param species_indexes
  * @param runs
  */
def writeAnnotatedRuns(index: Path, species_indexes: Path, runs: scala.List[AnnotatedRun]) = {
  index.toIO.asCsvWriter[AnnotatedRun](config.withHeader).write(runs)
  if (species_indexes != Path("/") && species_indexes.toIO.exists()) {
    val by_species = runs.groupBy(_.organism)
    for ((sp, rs) <- by_species) {
      Try {
        val p = species_indexes / (sp + ".tsv")
        p.toIO.asCsvWriter[AnnotatedRun](config.withHeader).write(rs)
        println(s"created per-species file for ${sp} at" + p.toIO.toScala.pathAsString)
      } match {
        case Failure(th) =>
          println(s"SPECIES FAILURE for ${sp}:")
          println(th)
        case _ =>
      }
    }
  }
}

/**
  * Builds an index of quantified GSE/GSM/SRR-s
 *
  * @param index index file name for samples
  * @param root root folder (/data/samples by default)
  * @return
  */
@main
def main(index: Path = Path("/data/samples/index.tsv"),
         root: Path = Path("/data/samples"),
         species_indexes: Path = Path("/"),
         rewrite: Boolean = false,
         ignoreFolders: Seq[Path] = Vector.empty
        ) = {
  val fl = index.toIO.toScala
  val rt = root.toIO.toScala
  val runs: List[AnnotatedRun] = rt.children.filter(_.isDirectory).flatMap {
    series =>
      println("============")
      println(s"SERIES = " +series.name)
      //(series.name, series.children.filter(_.isDirectory))
      val samples = series.children.filter(s => s.isDirectory && s.nonEmpty).toList
      samples.flatMap {
        case experiment if experiment.isDirectory && experiment.children.exists(f=> f.isDirectory &&
          f.children.exists(child=>child.name.contains("_transcripts_abundance.tsv"))
        ) =>
          println("Experiment = " + experiment.name)
          if(rewrite || !(experiment.children.exists(_.name == experiment.name + ".json") && experiment.children.exists(_.name == experiment.name + "_runs.tsv"))  )
          {
            println("cannot find " + experiment.name + ".json, creating it from scratch!")
            Try {
              val uname = experiment.name.toUpperCase
              if(uname.startsWith("PRJN") || uname.startsWith("SRX") || uname.startsWith("ERX"))
                geo.cli.MainCommand.fetchBioProject(experiment.name, key, experiment.pathAsString + "/" + experiment.name + ".json", experiment.pathAsString + "/" + experiment.name + "_runs.tsv", true)
              else geo.cli.MainCommand.fetchGSM(experiment.name, key, experiment.pathAsString + "/" + experiment.name + ".json", experiment.pathAsString + "/" + experiment.name + "_runs.tsv", true)

            } match {
              case Failure(th) => println("could not create Experiment because of: " + th.toString)
               case _ => println(experiment.name + ".json" + " successfully created!")
            }
          }
          //val sample = extractGSM(series, experiment)

          experiment.children.collect {
            case run if run.isDirectory && run.children.exists(_.name.contains("_transcripts_abundance.tsv")) =>
              println(series.name + "/" + experiment.name + "/" + run.name)
              FetchGEO(key).getAnnotatedRun(run.name, series.name,
                getPathIf(run)(_.name.contains("genes_abundance.tsv")),
                getPathIf(run)(_.name.contains("transcripts_abundance.tsv")),
                getPathIf(run)(_.name.contains("quant_")),
              )
          }
        case experiment =>
          println(s"Experiment ${experiment.name} does not seem to have SRR-s inside!")
          Nil
      }
  }.toList
  writeAnnotatedRuns(index, species_indexes, runs)
  println("INDEX SUCCESSFULLY CREATED at " + index.toIO.toScala.pathAsString)
}