#! /usr/local/bin/amm

import $exec.dependencies
import ammonite.ops._
import java.nio.file.Paths
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.


val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

def read_list(p: Path, header: Boolean = false): List[List[String]] = {
  p.toIO.unsafeReadCsv[List, List[String]](config.withHeader(header))
}

def merge(one: List[List[String]], two: List[List[String]]): List[List[String]] = {
  for {
    (a, b) <- one.zip(two)
  } yield if(a.head == b.head) a ++ b.tail else a ++ b
}

@main
def concat(where: Path, files: Path*): Path = {
  for(p <- files) read.lines(p).foreach(l=> write.append(where, l +"\n"))
  where
}

@main
def merge(where: Path, first: Path, other: Path*): Path = {
  val pathes = first :: other.toList
  pathes.map(p=>read_list(p)) match {
    case one::Nil =>
      where.toIO.writeCsv(one, config.withHeader(false))
      where
    case list =>
      val joined = list.reduce(merge)
      where.toIO.writeCsv(joined, config.withHeader(false))
      where
  }
}

def updateWithFolder(folder: Path, row: String): String = {
  (folder / row).toNIO.getFileName.toString
}

@main
def update_folder(from: Path, to: Path, folder: Path, indexes: Int*): Path = {
  val list = read_list(from)
  val where = indexes.toSet
  val updated: List[List[String]] = for{ row <- list  }
    yield
      for{ (r, i) <- row.zipWithIndex }
        yield if(where.contains(i) && r.trim!="" && r!="N/A") updateWithFolder(folder, r) else r
  to.toIO.writeCsv(updated, config.withHeader(false))
  to
}


@main
def copy_files_from_tsv(from: Path, destination: Path, indexes: Int*) = {
  val tsv = read_list(from, false)
  for{
    row <- tsv
    i <- indexes
  } cp.into(Path(row(i)), destination)
}

@main
def info() = {
  println("scripts to work with tsv")
}
