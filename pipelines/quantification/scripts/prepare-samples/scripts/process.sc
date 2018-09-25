#! /usr/local/bin/amm
import $exec.classes
import classes._

import ammonite.ops._
import java.nio.file.Paths
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.
import java.io.{File => JFile}

val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

def readIndexesJson(references: Path): Map[String, Indexes] = {
  val str: String = read(references)
  import io.circe.jackson.decode

  decode[Map[String, Indexes]](str) match {
    case Left(error) =>  throw error
    case Right(value) => value.map{ case (k, v) => k.toLowerCase -> v}
  }
}

@main
def main(samples: Path, references: Path, cache: Path): (Path, Path, Path) = process(samples, references, cache)

def read_list(p: Path, header: Boolean = false): List[List[String]] = {
  p.toIO.unsafeReadCsv[List, List[String]](config.withHeader(header))
}


@main
def write_sample(input_row: Path, sample_row: Path, to:Path, samples_folder: Path, gse: String, withHeaders: Boolean = false): Path = {
  val input = read_list(input_row)
  val gsm::tp::other = read_list(sample_row, false).head
  val toWrite: List[String] = input.head ++ (tp::FullSample.extractFiles(samples_folder, gse, gsm, other))
  if(withHeaders)
    to.toIO.writeCsv(List(toWrite), config.withHeader(FullSample.headers:_*))
  else to.toIO.writeCsv(List(toWrite), config.withHeader(false))
  to
}


@main
def process(samples: Path, references: Path, cache: Path): (Path, Path, Path) = {
  val indexes = readIndexesJson(references)
  processTSV(samples, cache, indexes)
}


def processTSV(samples: Path, cache_folder: Path, indexes: Map[String, Indexes]): (Path, Path, Path) = {
  val unsafeTSV = samples.toIO.unsafeReadCsv[Seq, Sample](config)
  val (valid_samples: Seq[Sample], invalid: Seq[Sample])= unsafeTSV.partition(s=>s.canQuantify(indexes))
  val valid = valid_samples.map(s=>s.toExtendedSample(indexes(s.species)))
  val dir = pwd //root / "data"

  val (cached, novel) =  valid.partition(p=>p.isInside(cache_folder, "sample.tsv"))
  val (cached_tsv: Path, novel_tsv: Path, invalid_tsv: Path) = (dir / "cached.tsv", dir / "novel.tsv", dir / "invalid.tsv")
  cached_tsv.toIO.writeCsv[ExtendedSample](cached, config.withHeader(false))
  novel_tsv.toIO.writeCsv[ExtendedSample](novel, config.withHeader(false))
  invalid_tsv.toIO.writeCsv[Sample](invalid, config.withHeader(false))
  (cached_tsv, novel_tsv, invalid_tsv)
}

@main
def update_from_json(path: Path, json: Path, indexes: (Int, String)*) = {
  val list = read_list(path)
  import JsonReader._
  val js = JsonReader.unsafeReadJson(json)
  val ins = indexes.map(_._1).toSet
  for {
    row <- list
  } yield row.zipWithIndex.map{
    case (r, i) if ins.contains(i) => js.unsafeReadString(r)
  }
}

@main
def update_from_json_column(from: Path, to: String, json_index: Int, indexes: (Int, String)*) = {
  val list = read_list(from)
  import JsonReader._
  val mp = indexes.toMap
  val ins = mp.keySet
  val lines = for {
    row <- list
    js = JsonReader.unsafeReadJson(Path(row(json_index)))
  }  yield row.zipWithIndex.map{ case (r, i) => if(ins.contains(i)) js.unsafeReadString(mp(i)) else r}
  new JFile(to).writeCsv(lines, config.withHeader(false))
}

def update_from_json_file(from: Path, to: String, json: Path, indexes: (Int, String)*) = {
  val list = read_list(from)
  import JsonReader._
  val js = JsonReader.unsafeReadJson(json)
  val mp = indexes.toMap
  val ins = mp.keySet
  val lines = for {
    row <- list
  } yield row.zipWithIndex.map{ case (r, i) => if(ins.contains(i))js.unsafeReadString(mp(i)) else r}
  new JFile(to).writeCsv(lines, config.withHeader(false))
}


@main
def info() = {
  println("prepare-samples script")
}
