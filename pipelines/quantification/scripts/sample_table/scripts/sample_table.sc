#! /usr/local/bin/amm
import $exec.dependencies
import geo.extras.SampleSummarizer
import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._
import io.circe.optics.JsonPath.root
import io.circe._
import io.circe.parser._
import geo.fetch._

import scala.util._
import better.files.File



lazy val key ="0a1d74f32382b8a154acacc3a024bdce3709"
implicit val ss = SampleSummarizer(FetchGEO(key))
import ss._

/**
 * Write indexes with runs to the disck
 * @param index
 * @param species_indexes
 * @param runs
 */
def writeAnnotatedRuns(index: Path, species_indexes: Path, runs: scala.List[AnnotatedRun], mergeTPMs: Boolean) = {
  println(s"WRITING RUNS TO ${index.toIO.toScala.pathAsString}")
  //pprint.pprintln(runs)

  index.toIO.asCsvWriter[AnnotatedRun](config.withHeader).write(runs).close()
  if (species_indexes != Path("/") && species_indexes.toIO.exists()) {
    val by_species = runs.groupBy(_.organism)
    for ((sp, rs) <- by_species) {
      Try {
        val p = species_indexes / (sp + ".tsv")
        p.toIO.asCsvWriter[AnnotatedRun](config.withHeader).write(rs).close()
        println(s"created per-species file for ${sp} at" + p.toIO.toScala.pathAsString)
      } match {
        case Failure(th) =>
          println(s"SPECIES FAILURE for ${sp}:")
          println(th)
        case _ =>
      }

      if(mergeTPMs){
        val trs = (species_indexes / (sp + "_transcripts.tsv")).toIO.toScala
        println(s"writing merged transcripts to ${trs.pathAsString}")
        mergeExpressions(trs, runs.map(r=>r.run ->File(r.transcripts)), "transcripts")
        val gs = (species_indexes / (sp + "_genes.tsv")).toIO.toScala
        println(s"writing merged gene expressions to ${trs.pathAsString}")
        mergeExpressions(gs, runs.map(r=>r.run ->File(r.genes)), "genes")
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
         ignore: Seq[Path] = Vector.empty
        ) = {
  implicit val f = FetchGEO(key)
  val fl = index.toIO.toScala
  val rt = root.toIO.toScala
  val ignoreFolders = ignore.map(_.toIO.toScala).toSet
  val rootChildren = rt.children.filter(c=> !ignoreFolders.contains(c))
  val runs: List[AnnotatedRun] = rt.children.filter(_.isDirectory).flatMap {
    series =>
      println("============")
      println(s"SERIES = " +series.name)
      //(series.name, series.children.filter(_.isDirectory))
      val samples = series.children.filter(s => s.isDirectory && s.nonEmpty && !ignoreFolders.contains(s)).toList
      samples.flatMap {
        case experiment if experiment.isDirectory && experiment.children.exists(f=> f.isDirectory &&
          f.children.exists(child=>child.name.contains("_transcripts_abundance.tsv"))
        ) =>
          processExperiment(series, experiment, rewrite)

        case experiment =>
          println(s"Experiment ${experiment.name} does not seem to have SRR-s inside!")
          Nil
      }
  }.toList
  writeAnnotatedRuns(index, species_indexes, runs, mergeTPMs = species_indexes.toIO.exists())
  println("INDEX SUCCESSFULLY CREATED at " + index.toIO.toScala.pathAsString)
}