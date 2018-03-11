#! /usr/local/bin/amm

import $exec.classes
import classes._
import $exec.tsv
import tsv._
import ammonite.ops._
import java.nio.file.Paths

import io.circe.Json
import kantan.csv._
import kantan.csv.ops._     // Enriches types with useful methods.


val result = FullSample.getLibRatio(Path("/pipelines/test/lib_format_counts.json"))
//val js = "/pipelines/test/lib_format_counts.json"
pprint.pprintln(result)

object Extractor extends JSONReader

@main
def update_files(path: Path, json: Path, indexes: (Int, String)*) = {
  val list = read_list(path)
  val js = Extractor.unsafeReadJson(json)
  val ins = indexes.map(_._1).toSet
  for {
    row <- list
  } yield row.zipWithIndex.map{
    case (r, i) if ins.contains(i) => Extractor.unsafeReadField(js, r)
  }
}