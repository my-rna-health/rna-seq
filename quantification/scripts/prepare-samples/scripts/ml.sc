#! /usr/local/bin/amm
import $exec.classes
import classes._

import ammonite.ops._
import java.nio.file.Paths
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.
import java.io.{File => JFile}


implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

@main
def summarize(p: Path, destination: Path, header: Boolean = false) = {
  val groups = FullSample.read_samples(p)(config.withHeader(header)).groupBy(s=>s.sample_group)
  for{
    (title, group) <- groups
    }
  {
    val d = destination / (title + ".tsv")
    val headers: List[String] = group.head.getExpressions(config.withHeader(true)).map(e=>e.Name)
    val len = headers.size
    val rows = group.map{ g=>
      val exp = g.getExpressions
      require(exp.length == len, s"for the same salmon index there should be same number of transcripts, currently we have ${len} vs ${exp.size}")
      g.toList ++ g.getExpressions.map(e=>e.TPM.toString)
    }
    d.toIO.writeCsv(rows, config.withHeader(FullSample.headers ++ headers:_*))
  }
}

@main
def info() = {
  println("prepare machine learning")
}
