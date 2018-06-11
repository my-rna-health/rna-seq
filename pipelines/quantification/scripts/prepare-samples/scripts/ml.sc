#! /usr/local/bin/amm
import geny.Generator
import $exec.classes
import classes._

import scala.collection.immutable
import ammonite.ops._
import java.nio.file.Paths
import java.io.{File => JFile}
import better.files._
import File._


import ammonite.main.Router.main

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

protected def write_expression_row(where: Path, quant: Seq[(String, Double)], sample: String) = {
  /*
  val (headers, values) = quant.foldLeft(("sample", sample)){
    case ((h, v), (name, value)) => ( h + "\t" + name, v + "\t" + value )
  }
  write.append(where, headers + "\n")
  write.append(where, values + "\n")
  */
  write.append(where, "samples\\transcripts")
  for((t, _) <- quant) write.append(where, s"\t ${t}")
  write.append(where, s"\n${sample}")
  for((_, v) <- quant) write.append(where, s"\t ${v}")
  write.append(where, s"\n")
  where
}

protected def write_expression_column(where: Path, quant: Seq[(String, Double)], sample: String) = {
  write.append(where, s"transcripts\\samples\t${sample}\n")
  for((name, value) <- quant){
    write.append(where, s"${name}\t${value}\n")
  }
  where
}

@main
def write_sample_expressions(gsm: Path): Option[(Path, Path)] = {
  val row = gsm / "expression_row.tsv"
  val column = gsm / "expression_column.tsv"
  val q = sampleQuantPath(gsm)
  if(! (exists! q)) {
    println(s"quantification for ${gsm} does not exist!")
    None
  } else {
      val name = gsm.name
      val quant = SalmonExpressions.read_quants(q).map(q=>q.Name -> q.TPM)//.read_named_TPMs(q)
      println(s"read ${name}'s expressions")
      if(exists! row) println(s"expression row is already written to ${row} !") else {
        write_expression_row(row, quant, gsm.name)
        println(s"expression rows are written for ${name}")
      }
      if(exists! column) println(s"expression row is already written to ${column} !") else {
        write_expression_column(column, quant, name)
        println(s"expressions columns are written for ${name}")
      }
    }
    Some((row, column))
}

protected def writeRows(where: Path, rows: Seq[Path]) = {
  write.write(where, "")
  if(rows.isEmpty){
    println("no quantified samples inside series!")
    None
  } else{
   val header = (read.lines.iter! rows.head).head + "\n"
   write.append(where, header)
   for{r <- rows
     str <- (read.lines.iter! r).drop(1)
   } write.append(where, str + "\n")
   println(s"rows written to ${where}")
  }
  Some(where)
}

protected def writeCols(where: Path, columns: Seq[Path]) = {
  write.write(where, "")
  if(columns.isEmpty){
    println("no quantified samples inside series!")
    None
  } else {
    val first = columns.head
    val transLines = (read.lines.iter! first).map(l=>l.takeWhile(c=>c!='\t'))
    val iters: Seq[Iterator[String]] = columns.map(c=>File(c.toString).lineIterator)
    for(l <- transLines){
     write.append(where, l)
     for(i <- iters) {
           write.append(where, i.next().dropWhile(c=>c!='\t'))
     }
     write.append(where, "\n")
    }
    println(s"columns written to ${where}")
  }
  Some(where)
}

@main
def write_series_expressions(gse: Path) = {
  val gsms: LsSeq = ls! gse
  val name = gse.name
  val (samples, unfinished) = gsms.filter(p=>p.isDir).partition(g=> exists! sampleQuantPath(g))
  for(u <- unfinished) println(s"cannot find quantification for ${u}")
  val quants: Seq[(Path, Path)] = samples.map { s => write_sample_expressions(s) }.collect{ case Some( (row, column) ) => (row, column)}
  println("---------------------------------------")
  println(s"writing rows and columns for ${name}")
  val rows = gse / s"expressions_rows_${name}.tsv"
  writeRows(rows, quants.map(_._1))
  val cols = gse / s"expressions_columns_${name}.tsv"
  writeCols(cols, quants.map(_._2))
  (rows, cols)
}


@main
def info(): Unit = {
  println("scripts to prepare expression data for machine learning")
}