#! /usr/local/bin/amm
import $exec.classes
import classes._

import ammonite.ops._
import java.nio.file.Paths
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.


val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

def readIndexesJson(references: Path): Map[String, Indexes] = {
  val str: String = read(references)
  import io.circe.jackson.decode

  decode[Map[String, Indexes]](str) match {
    case Left(error) =>  throw error
    case Right(value) => value
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
  val (valid_samples: Seq[Sample], invalid: Seq[Sample])= samples.toIO.unsafeReadCsv[Seq, Sample](config).partition(s=>s.canQuantify(indexes))
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
def info() = {
  println("prepare-samples script")
}
