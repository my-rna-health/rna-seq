#! /usr/local/bin/amm
import $exec.classes
import classes._

import ammonite.ops._
import java.nio.file.Paths
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.
import java.io.{File => JFile}


val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)


def read_samples(p: Path, header: Boolean = false): List[FullSample] = {
  val rows: List[List[String]] = p.toIO.unsafeReadCsv[List, List[String]](config.withHeader(header))
  rows.map(row => FullSample.fromList(row))
}

def read_quants(p: Path, header: Boolean = false): List[SalmonExpressions] = {
  p.toIO.unsafeReadCsv[List, SalmonExpressions](config.withHeader(header))
}


//val p = Path("/pipelines/batches/novel.tsv")
//val q = Path("/pipelines/test/quant.sf")
//pprint.pprintln(read_samples(p))
//pprint.pprintln(read_quants(q).take(10))


