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

def read_vector(p: Path, header: Boolean = false): Vector[Vector[String]] = {
  p.toIO.unsafeReadCsv[Vector, Vector[String]](config.withHeader(header))
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


//TODO: rename
@main
def concat_one_header(where: Path, files: Path*): Path = {
  //we assume that second and later files have headers, so we skip first row there
  for(p <- files) read.lines(p).tail.foreach(l=> write.append(where, l +"\n"))
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
def copy_indexed_files(from: Path, destination: String, indexes: Int*) = {
  val tsv = read_list(from, false)
  val d = Path(destination)
  if(!exists(d)) mkdir(d)
  for{
    row <- tsv
    i <- indexes
  } {
    cp.into(Path(row(i)), d)
  }
}

@main
def transpose(from: Path, to: Path) = {
  val vec = read_vector(from).transpose[String]
  to.toIO.writeCsv(vec, config.withoutHeader)
  to
}

/*
def prefixesOrSubpathes(p: Path, prefix: String *) = {
  val subs: Seq[String] = prefix.takeWhile(p=>p.endsWith("/"))
  //val prefixes = prefix.sk
}
*/

@main
def copy_prefixed_files(from: Path, destination: String, prefix_index_1: Int, prefix_index_2: Int, prefix_index_3: Int,  indexes: Int*) = {
  //example: docker run -v /pipelines:/pipelines quay.io/comp-bio-aging/prepare-samples tsv.sc copy_prefixed_files /pipelines/samples/batches/cross-species-1/novel.tsv /pipelines/samples/bes/cross-species-1/quant 2 8 0 23 24
  val tsv = read_list(from, false)
  val d = Path(destination)
  if(!exists(d)) mkdir(d)
  for{
    row <- tsv
    prefix_1 = if(prefix_index_1 >= 0) row(prefix_index_1) + "_" else ""
    prefix_2 = if(prefix_index_2 >= 0) row(prefix_index_2) + "_" else "_"
    prefix_3 = if(prefix_index_3 >= 0) row(prefix_index_3) + "_" else "_"
    i <- indexes
  } {
    val p = Path(row(i))
    val name = s"${prefix_1}${prefix_2}${prefix_3}".replace(" ", "_").replace("__", "_") + p.segments.last
    cp.over(p, d / name)
  }
}


@main
def info() = {
  println("scripts to work with tsv")
}
