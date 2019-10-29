#! /usr/local/bin/amm
import $exec.dependencies
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
import geo.cli._

import scala.util._
import better.files.File
import kantan.csv.{CsvConfiguration, rfc}
implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)


def fix(proj: better.files.File, root: better.files.File, key: String = "0a1d74f32382b8a154acacc3a024bdce3709") = {
  val f = FetchGEO(key)
  val name = proj.name
  println(s"fix for $name")
  val place = proj.moveToDirectory(root)
  val children = place.children.filter(c =>  c.isDirectory).toVector
  val bio =f.getBioProject(name)
  val srx = bio.experimentIds
  for(s <- srx){
    val (e, runs) = f.runsFromExperiment(f.getExperiment(s))
    val run_ids = runs.map(_.run.Run).toSet
    val rs = children.filter(c=>run_ids.contains(c.name))
    if(rs.nonEmpty){
      println(s"${name}: fix for SRX $s")
      val sf = (place / s).createDirectory()
      rs.foreach(_.moveToDirectory(sf))
      val runsPath = (sf / (sf.name + "_runs.tsv")).pathAsString
      val o = (sf / (sf.name + ".json")).pathAsString
      CommandsBioProject.fetchExperiment(s, key, o, runsPath, true)
    }
    //val rs = runs.map(r => Run(proj, s, r.run.Run, r.sample.TaxID, r.run.Experiment, r.sample.Model, r.library.LibraryStrategy, r.library.LibrarySelection, r.library.LibrarySource, "", e.experiment.))
    //(sf / (s + "_runs.tsv")).toJava.asCsvWriter[Run](config.withHeader).write()
  }
}

@main def fix_folders(root: Path = Path("/data/samples/species"), key: String = "0a1d74f32382b8a154acacc3a024bdce3709") = {
  val projects = root.toIO.toScala.children.filter(_.name.startsWith("PRJN")).toList
  println("starting the fix for wrong folders inside projects")
  val f = FetchGEO(key)
  for(proj <- projects){
    val name = proj.name
    val bio =f.getBioProject(name)
    val srx = bio.experimentIds
    val children = proj.children.toList
    val wrong = children.filter(_.name.contains("SRR")).toList
    //val run_ids = wrong.map(_.name)
    for(s <- srx){
      val (e, runs) = f.runsFromExperiment(f.getExperiment(s))
      val child_run_ids = runs.map(_.run.Run).toSet
      val rs = children.filter(c=>child_run_ids.contains(c.name)).toList
      if(rs.nonEmpty){
        println(s"${name}: fix for SRX $s with runs: [${rs.map(_.name).mkString(",")}]")
        val sf = (proj / s).createDirectoryIfNotExists()
        rs.foreach(_.moveToDirectory(sf))
        val runsPath = (sf / (sf.name + "_runs.tsv")).pathAsString
        val o = (sf / (sf.name + ".json")).pathAsString
        CommandsBioProject.fetchExperiment(s, key, o, runsPath, true)
      }
      //val rs = runs.map(r => Run(proj, s, r.run.Run, r.sample.TaxID, r.run.Experiment, r.sample.Model, r.library.LibraryStrategy, r.library.LibrarySelection, r.library.LibrarySource, "", e.experiment.))
      //(sf / (s + "_runs.tsv")).toJava.asCsvWriter[Run](config.withHeader).write()
    }
  }
}

@main def bioproject_fix(bioprojects: Path = Path("/data/samples/species/bioprojects"), root:  Path = Path("/data/samples/species/")) = {
  val r  = root.toIO.toScala
  val prs = bioprojects.toIO.toScala.children.filter(_.isDirectory).toVector
  println(s"FIXING ISSUES IN $prs")
  for(p <- prs) fix (p, r)
}
