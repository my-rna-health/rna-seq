import scala.collection.immutable
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


protected def writeColumns(summary: Path, expressions: Seq[(String, Seq[SalmonExpressions])]): Path = {
  val transcripts = expressions.head._2.map(e=>e.Name).zipWithIndex
  val sampleNames = expressions.map(_._1).toList
  val header: String = ("transcripts/samples"::sampleNames).mkString("\t").dropRight(1) + "\n"
  write.append(summary, header)
  for{
    (name, i) <- transcripts
  }{
    val str = name + "\t" + expressions.map(_._2(i).TPM).mkString("\t").dropRight(1) + "\n"
    write.append(summary, str)
  }
  summary
}

protected def writeRows(summary: Path, expressions: Seq[(String, Seq[SalmonExpressions])]): Path = {
  val transcriptNames = expressions.head._2.map(e=>e.Name).toList
  val headers = transcriptNames.foldLeft("samples/transcripts"){ case (acc, el) => acc + "\t" + el} + "\n"
  write.append(summary, headers)
  for((s, e)<- expressions)
    {
      val str = e.foldLeft(s){case (acc, el) => acc + "\t" + el.TPM} + "\n"
      write.append(summary, str)
    }
  summary
}

/*
protected def getSeriesExpressions(p: Path): Seq[(String, Seq[SalmonExpressions])]= {
  val gsms: LsSeq = ls! p
  val (samples, unfinished) = gsms.partition(p=> exists! sampleQuantPath(p))
  for(p <- unfinished) println(s"cannot find quantification for ${p}")
  samples.toVector.map(p => (p.name, SalmonExpressions.read_quants(sampleQuantPath(p))(config)))
}
*/

protected def getSeriesExpressions(p: Path): Vector[(String, Vector[(String, Double)])] = {
  val gsms: LsSeq = ls! p
  val (samples, unfinished) = gsms.partition(p=> exists! sampleQuantPath(p))
  for(p <- unfinished) println(s"cannot find quantification for ${p}")
  samples.toVector.map { p =>
    val qPath = sampleQuantPath(p)
    (p.name, SalmonExpressions.read_named_TPMs(qPath)(config))
  }
}

@main
def series_expressions(p: Path, summary: Path, columns: Boolean): Path = {
  val expressions = getSeriesExpressions(p)
  println("TO THE END")
  p
  //if(columns) writeColumns(summary, expressions) else writeRows(summary, expressions)
}

@main
def info(): Unit = {
  println("scripts to prepare expression data for machine learning")
}


//val p = Path("/pipelines/RESULTS/GSE75192/")

//series_expressions(p, p / "columns.tsv", true)
//series_expressions(p, p / "rows.tsv", false)