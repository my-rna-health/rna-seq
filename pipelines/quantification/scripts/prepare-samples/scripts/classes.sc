import $exec.dependencies
import ammonite.ops._
import java.io.{File => JFile}
import java.nio.file.{Paths, Path => JPath}
//import scala.collection.immutable._

import io.circe.{Decoder, Json}
import io.circe.generic.JsonCodec
import kantan.csv._         // All kantan.csv types.
import kantan.csv.ops._     // Enriches types with useful methods.


@JsonCodec case class Indexes(salmon: String, transcriptome: String, gtf: String)

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

  def filename(str: String) = Path(str).segments.last

}

object ExtendedSample {
  lazy val headers = List("GSM",	"GSE",	"Species",	"Sequencer",	"Type",
    "Sex",	"Age",	"Tissue",	"Extracted molecule",
    "Strain",	"Comments", "salmon", "transcriptome", "gtf")

  implicit val extendedSampleCodec: HeaderCodec[ExtendedSample] = HeaderCodec.caseCodec(
    "GSM",	"GSE",	"Species",	"Sequencer",	"Type",
    "Sex",	"Age",	"Tissue",	"Extracted molecule",
    "Strain",	"Comments", "salmon", "transcriptome", "gtf")(ExtendedSample.apply)(ExtendedSample.unapply)
}

case class ExtendedSample(gsm: String,	gse: String,	species: String,
                          sequencer: String, sample_type: String,	sex: String,
                          age: String,	tissue: String,	extracted_molecule: String,
                          strain: String,	comments: String, salmon: String, transcriptome: String, gtf: String //gtf is optional
                         ) extends PipelineSample

object Sample {

  implicit val sampleCodec: HeaderCodec[Sample] = HeaderCodec.caseCodec(
    "GSM",	"GSE",	"Species",	"Sequencer",	"Type",
    "Sex",	"Age",	"Tissue",	"Extracted molecule",
    "Strain",	"Comments")(Sample.apply)(Sample.unapply)
}


case class Sample(gsm: String,	gse: String,	species: String,
                  sequencer: String, sample_type: String,	sex: String,
                  age: String,	tissue: String,	extracted_molecule: String,
                  strain: String,	comments: String) extends PipelineSample
{


  def toExtendedSample(index: Indexes): ExtendedSample = {
    ExtendedSample(gsm, gse, species, sequencer, sample_type, sex, age, tissue, extracted_molecule, strain, comments,
      index.salmon, index.transcriptome, index.gtf)
  }

}

/*
A = automatic

I = inward
O = outward
M = matching
+
S = stranded
U = unstranded
+
F
R
 */
//def libType: String = "A"

object JsonReader extends JsonReader
{
}

trait JsonReader {

  def unsafeReadJson(references: Path): Json = {
    import io.circe.jackson.decode
    import io.circe.jackson.parse
    val str: String = read(references)
    val json = parse(str).toOption.get
    json
  }

  implicit class JsonExt(json: Json)
  {
    def unsafeReadField[A](field: String)(implicit d: Decoder[A]) = json.hcursor.get[A](field).toOption.get

    def unsafeReadString(field: String)(implicit dstring: Decoder[String], dnum: Decoder[Double]) = {
      json.hcursor.get[String](field)(dstring).getOrElse(json.hcursor.get[Double](field)(dnum).map(d=>d.toString).toOption.get)
    }

  }

}

trait FolderExtractor {

  def extractNames(list: List[String]): List[String] = list.map(s=> Paths.get(s).getFileName.toString)
  def toFolder(folder: Path, value: String): String = if(value.trim=="" || value =="N/A") "" else (folder / value).toString()
  def toFolder(folder: Path, list: List[String]): List[String] = list.map(s=>toFolder(folder, s))

}


case class SalmonExpressions(Name: String,	Length: Int,	EffectiveLength: Double,	TPM: Double,	NumReads: Double)

object SalmonExpressions {
  val headers = Seq("Name",	"Length",	"EffectiveLength",	"TPM",	"NumReads")
  implicit val salmonExpressionsCodec: HeaderCodec[SalmonExpressions] = HeaderCodec.caseCodec("Name",	"Length",	"EffectiveLength",	"TPM",	"NumReads")(SalmonExpressions.apply)(SalmonExpressions.unapply)

  def read_quants(p: Path)(implicit config: CsvConfiguration): Seq[SalmonExpressions] = {
    p.toIO.unsafeReadCsv[Vector, SalmonExpressions](config)
  }

}

trait SalmonSample {
  def quant: String
  def expressions: String
  def compatible_fragment_ratio: String
  def expected_format: String

}

object FullSample extends FolderExtractor with JsonReader {

  def fromList(list: List[String]) = list match {
    case  gsm::	gse::	species:: sequencer:: sample_type::	sex::  age::	tissue::	extracted_molecule::  strain::	comments:: //Sample
      salmon:: transcriptome:: gtf:: //ExtendedSample
      reads_type::forward_read_raw:: reverse_read_raw:: sra::  forward_read_cleaned:: reverse_read_cleaned:: quality_html:: quality_json::
      quant:: expressions:: lib_format:: expected_format:: compatible_fragment_ratio :: _ =>
        FullSample(gsm,	gse,	species, sequencer, sample_type,	sex,  age,	tissue,	extracted_molecule,  strain,	comments, //Sample
          salmon, transcriptome, gtf, //ExtendedSample
          reads_type, forward_read_raw, reverse_read_raw, sra,  forward_read_cleaned, reverse_read_cleaned, quality_html, quality_json,
          quant, expressions, lib_format, expected_format, compatible_fragment_ratio)
    case _ => throw new Exception(s"list does not have enough fields! ${list.mkString(" ")}")
  }

  def read_samples(p: Path)(implicit config: CsvConfiguration): List[FullSample] = {
    val rows: List[List[String]] = p.toIO.unsafeReadCsv[List, List[String]](config)
    rows.map(row => FullSample.fromList(row))
  }


  def getLibRatio(path: Path): (String, Double) = {
    val json = unsafeReadJson(path)
    println()
    val ratio = json.unsafeReadField[Double]("compatible_fragment_ratio")
    println(ratio)
    val lib = json.unsafeReadField[String]("expected_format")
    println(lib)
    (lib, ratio)
  }

  /*
  lazy val sampleFilesHeaders = List(
    "forward_read_raw", "reverse_read_raw", "sra", //basic
    "forward_read_cleaned", "reverse_read_cleaned", "quality_html", "quality_json", //fastp
    "quant", "expressions",  "lib_format",  "expected_format", "compatible_fragment_ratio" //quantification
  )
  */

  //val headers =  ExtendedSample.headers ++sampleFilesHeaders

  val headers = Seq("GSM",	"GSE",	"Species",
    "Sequencer",	"Type",  "Sex",
    "Age",	"Tissue",	"Extracted molecule",
    "Strain",	"Comments", //Sample
    "salmon", "transcriptome", "gtf", //indexes
    "reads_type","forward_read_raw", "reverse_read_raw", "sra", //basic
    "forward_read_cleaned", "reverse_read_cleaned", "quality_html", "quality_json", //fastp
    "quant", "expressions",  "lib_format",  "expected_format", "compatible_fragment_ratio" //quantification)
  )

  def extractFiles(samples_folder: Path, gsm: String, gse: String, files: List[String]) = {
    val destination = samples_folder / gse / gsm
    val names: List[String] = extractNames(files)
    val reads = names.head::names.tail.head::Nil
    val toWrite: List[String] = toFolder(destination / "raw", names) ++
      toFolder(destination / "cleaned", reads.map(_.replace(".fastq", "_cleaned.fastq"))) ++
      toFolder(destination / "report" , List("fastp.html", "fastp.json")) ++
      (toFolder(destination, "transcripts_quant")::
        toFolder(destination / "transcripts_quant" , List("quant.sf", "lib_format_counts.json", "lib_format_counts.json" ,"lib_format_counts.json")))

    toWrite
  }
  //implicit val fullSampleDecoder: RowDecoder[FullSample] = HeaderDecoder.defaultHeaderDecoder[FullSample].fromHeader(FullSample.headers).toOption.get

}


case class FullSample(gsm: String,	gse: String,	species: String,
                      sequencer: String, sample_type: String,	sex: String,
                      age: String,	tissue: String,	extracted_molecule: String,
                      strain: String,	comments: String, //Sample
                      salmon: String, transcriptome: String, gtf: String, //ExtendedSample
                      reads_type: String,  forward_read_raw: String, reverse_read_raw: String, sra: String,
                      forward_read_cleaned: String, reverse_read_cleaned: String, quality_html: String, quality_json: String,
                      quant: String, expressions: String, lib_format: String, expected_format: String, compatible_fragment_ratio: String
                     ) extends PipelineSample with SalmonSample {

  def getExpressions(implicit config: CsvConfiguration): Seq[SalmonExpressions] = SalmonExpressions.read_quants(Path(expressions))(config)


  def sample_group =s"${species}-${filename(transcriptome)}-${extracted_molecule}".replace(" ", "_")

  lazy val toList: List[String] = gsm::	gse::	species:: sequencer:: sample_type::	sex::  age::	tissue::	extracted_molecule::  strain::	comments:: //Sample
    salmon:: transcriptome:: gtf:: //ExtendedSample
    forward_read_raw:: reverse_read_raw:: sra::  forward_read_cleaned:: reverse_read_cleaned:: quality_html:: quality_json::
    quant:: expressions:: lib_format:: expected_format:: compatible_fragment_ratio ::Nil

}