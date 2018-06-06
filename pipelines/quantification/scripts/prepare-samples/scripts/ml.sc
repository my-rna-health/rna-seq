#! /usr/local/bin/amm
import $exec.classes
import classes._
import scala.collection.immutable
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


protected def writeRows(where: Path, samples: Seq[(String, Path)], log: Boolean): Path = if(samples.isEmpty) {
    println(s"no samples discovered inside of the series, ${where} will be empty!")
    where
  } else {
  val headers: String = SalmonExpressions.read_names(samples.head._2).foldLeft("samples/transcripts"){ case (acc, el) => acc + "\t" + el} + "\n"
  write.append(where, headers)
  if(log) println(s"headers are written to ${where}")
  for{
    (s, file) <- samples
  }{
    val e = SalmonExpressions.read_named_TPMs(file)
    val str = e.foldLeft(s){case (acc, el) => acc + "\t" + el._2} + "\n"
    write.append(where, str)
    if(log) println(s"expressions of ${s} from ${file} are written to ${where}")
  }
  where
}

protected def getSeriesExpressions(p: Path): Seq[(String, Path)] = {
  val gsms: LsSeq = ls! p
  val (samples, unfinished) = gsms.partition(p=> exists! sampleQuantPath(p))
  for(p <- unfinished) println(s"cannot find quantification for ${p}")
  samples.map { p =>
    (p.name, sampleQuantPath(p))
  }
}

@main
def series_expressions(p: Path, log: Boolean = false): Path = {
  val exp = p / (s"expressions_rows_${p.name}.tsv")
  val expressions = getSeriesExpressions(p)
  if(log) println(s"${expressions.size} samples found in ${p}")
  writeRows(exp, expressions, log)
  println(s"expression rows of ${p.name} written to ${exp}")
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