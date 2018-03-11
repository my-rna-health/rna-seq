#! /usr/local/bin/amm

import $exec.classes
import classes._
import classes.FullSample
import $exec.tsv
import tsv._
import ammonite.ops._
import java.nio.file.Paths
import java.io.{File => JFile}

import io.circe.Json
import kantan.csv._
import kantan.csv.ops._     // Enriches types with useful methods.


def read_samples(p: Path, header: Boolean = false): List[FullSample] = {
  p.toIO.unsafeReadCsv[List, FullSample](config.withHeader(header))
}

val p = Path("/pipelines/batches/novel.tsv")
val samples = read_samples(p, false)
pprint.pprintln(samples)