#! /usr/local/bin/amm
import $exec.classes
import classes._

import ammonite.ops._
import java.nio.file.Paths
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.
import java.io.{File => JFile}

implicit val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

protected def sampleQuantPath(p: Path): Path = p / "transcripts_quant" /  "quant.sf"

@main
def summarize(p: Path, destination: Path, header: Boolean = false) = {
  val groups = FullSample.read_samples(p)(config.withHeader(header)).groupBy(s=>s.sample_group)
  for{
    (title, group) <- groups
    }
  {
    val d = destination / (title + ".tsv")
    val headers = group.head.getExpressions(config.withHeader(true)).map(e=>e.Name)
    val len = headers.size
    val rows: List[List[String]] = group.map{ g=>
      val exp = g.getExpressions
      require(exp.length == len, s"for the same salmon index there should be same number of transcripts, currently we have ${len} vs ${exp.size}")
      g.toList ++ g.getExpressions.map(e=>e.TPM.toString)
    }
    d.toIO.writeCsv(rows, config.withHeader(FullSample.headers ++ headers:_*))
  }
}


protected def getSeriesExpressions(p: Path): Seq[(String, Seq[SalmonExpressions])]= {
  val gsms: LsSeq = ls! p
  val (samples, unfinished) = gsms.partition(p=> exists! sampleQuantPath(p))
  for(p <- unfinished) println(s"cannot find quantification for ${p}")
  samples.toVector.map(p => (p.name, SalmonExpressions.read_quants(sampleQuantPath(p))(config)))
}

protected def writeColumns(path: Path, expressions: Seq[(String, Seq[SalmonExpressions])]): Path = {
  val transcripts = expressions.head._2.map(e=>e.Name).zipWithIndex
  val sampleNames = expressions.map(_._1).toList
  val header: String = ("transcripts/samples"::sampleNames).mkString("\t").dropRight(1) + "\n"
  write.append(path, header)
  for{
    (name, i) <- transcripts
  }{
    val str = name + "\t" + expressions.map(_._2(i).TPM).mkString("\t").dropRight(1) + "\n"
    write.append(path, str)
  }
  path
}

protected def writeRows(path: Path, expressions: Seq[(String, Seq[SalmonExpressions])]): Path = {
  val transcriptNames = expressions.head._2.map(e=>e.Name).toList
  val headers: String = ("samples/transcripts"::transcriptNames).mkString("\t").dropRight(1) + "\n"
  write.append(path, headers)
  for((s, e)<- expressions)
    {
      val str = s + "\t" + e.map(_.TPM).mkString("\t").dropRight(1) + "\n"
      write.append(path, str)
    }
  path
}

@main
def series_expressions(p: Path, columns: Boolean): Path = {
  val summary = p / "summary.tsv"
  val expressions = getSeriesExpressions(p)
  if(columns) writeColumns(p, expressions) else writeRows(p, expressions)
}

@main
def info(): Unit = {
  println("scripts to prepare expression data for machine learning")
}