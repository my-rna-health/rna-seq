#! /usr/local/bin/amm

import $exec.dependencies
import ammonite.ops._
import better.files._
import java.io.{File => JFile}
import io.circe.generic.JsonCodec, io.circe.syntax._
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.


trait PipelineSample {
  def gsm: String
  def gse: String
  def species: String
  def sequencer: String
  def sample_type: String
  def sex: String
  def age: String
  def tissue: String
  def extracted_molecule: String
  def strain: String
  def comments: String

  def isInside(folder: Path): Boolean = {
    val p = folder / gse / gsm
    exists(p) && p.toIO.isDirectory
  }

  def isInside(folder: Path, sub: String): Boolean = {
    val p = folder / gse / gsm / sub
    exists(p)
  }

  def canQuantify(indexes: Map[String, Indexes]): Boolean = {
    indexes.contains(species) && indexes(species).salmon != ""
  }

}

case class FullSample(gsm: String,	gse: String,	species: String,
                          sequencer: String, sample_type: String,	sex: String,
                          age: String,	tissue: String,	extracted_molecule: String,
                          strain: String,	comments: String, salmon: String, transcriptome: String, gtf: String //gtf is optional
                     ) extends PipelineSample

case class Sample(gsm: String,	gse: String,	species: String,
                  sequencer: String, sample_type: String,	sex: String,
                  age: String,	tissue: String,	extracted_molecule: String,
                  strain: String,	comments: String) extends PipelineSample
{


  def toFullSample(index: Indexes): FullSample = {
        FullSample(gsm, gse, species, sequencer, sample_type, sex, age, tissue, extracted_molecule, strain, comments,
           index.salmon, index.transcriptome, index.gtf)
  }

}

implicit val sampleCoder: HeaderCodec[Sample] = HeaderCodec.caseCodec(
  "GSM",	"GSE",	"Species",	"Sequencer",	"Type",
  "Sex",	"Age",	"Tissue",	"Extracted molecule",
  "Strain",	"Comments")(Sample.apply)(Sample.unapply)

implicit val fullSampleCoder: HeaderCodec[FullSample] = HeaderCodec.caseCodec(
  "GSM",	"GSE",	"Species",	"Sequencer",	"Type",
  "Sex",	"Age",	"Tissue",	"Extracted molecule",
  "Strain",	"Comments", "salmon", "transcriptome", "gtf")(FullSample.apply)(FullSample.unapply)

val config: CsvConfiguration = rfc.withCellSeparator('\t').withHeader(true)

import io.circe.syntax._
@JsonCodec case class Indexes(salmon: String, transcriptome: String, gtf: String)

def readJson(references: Path): Map[String, Indexes] = {
  val str: String = read(references)
  import io.circe.jackson.decode

  decode[Map[String, Indexes]](str) match {
    case Left(error) =>  throw error
    case Right(value) => value
  }
}

@main
def main(samples: Path, references: Path, cache: Path): (Path, Path, Path) = process(samples, references, cache)


def merge(one: List[List[String]], two: List[List[String]]): List[List[String]] = {
  for {
    (a, b) <- one.zip(two)
  } yield if(a.head == b.head) a ++ b.tail else a ++ b
}

def read_list(p: Path, header: Boolean = false): List[List[String]] = {
  p.toIO.unsafeReadCsv[List, List[String]](config.withHeader(header))
}


@main
def concat(where: Path, files: Path*) = {
  val f = where
  for(p <- files) read.lines(p).foreach(l=> write.append(f, l +"\n"))
  f
}

@main
def merge(where: Path, first: Path, other: Path*): Path = {
  val pathes = first :: other.toList
  val list: List[List[List[String]]] = pathes.map(p=>read_list(p))
  val joined = list.reduce(merge)
  where.toIO.writeCsv(joined, config.withHeader(false))
  where
}

def updateWithFolder(folder: Path, row: String): String = {
  (folder / Path(row).toNIO.getFileName.toString).toString
}

@main
def update_folder(from: Path, to: Path, folder: Path, indexes: Int*) = {
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
def process(samples: Path, references: Path, cache: Path): (Path, Path, Path) = {
  val indexes = readJson(references)
  processTSV(samples, cache, indexes)
}


def processTSV(samples: Path, cache_folder: Path, indexes: Map[String, Indexes]): (Path, Path, Path) = {
  val (valid_samples: Seq[Sample], invalid: Seq[Sample])= samples.toIO.unsafeReadCsv[Seq, Sample](config).partition(s=>s.canQuantify(indexes))
  val valid = valid_samples.map(s=>s.toFullSample(indexes(s.species)))
  val dir = pwd //root / "data"

  val (cached, novel) =  valid.partition(p=>p.isInside(cache_folder, "sample.tsv"))
  val (cached_tsv: Path, novel_tsv: Path, invalid_tsv: Path) = (dir / "cached.tsv", dir / "novel.tsv", dir / "invalid.tsv")
  cached_tsv.toIO.writeCsv[FullSample](cached, config.withHeader(false))
  novel_tsv.toIO.writeCsv[FullSample](novel, config.withHeader(false))
  invalid_tsv.toIO.writeCsv[Sample](invalid, config.withHeader(false))
  (cached_tsv, novel_tsv, invalid_tsv)
}