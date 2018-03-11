import $exec.dependencies
import ammonite.ops._
import java.io.{File => JFile}
import java.nio.file.{Paths, Path => JPath}
import io.circe.{Decoder, Json}
import io.circe.generic.JsonCodec
import kantan.csv._

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

}

object ExtendedSample {
  lazy val headers = List("GSM",	"GSE",	"Species",	"Sequencer",	"Type",
    "Sex",	"Age",	"Tissue",	"Extracted molecule",
    "Strain",	"Comments", "salmon", "transcriptome", "gtf")
}

case class ExtendedSample(gsm: String,	gse: String,	species: String,
                          sequencer: String, sample_type: String,	sex: String,
                          age: String,	tissue: String,	extracted_molecule: String,
                          strain: String,	comments: String, salmon: String, transcriptome: String, gtf: String //gtf is optional
                         ) extends PipelineSample

implicit val extendedSampleCodec: HeaderCodec[ExtendedSample] = HeaderCodec.caseCodec(
  "GSM",	"GSE",	"Species",	"Sequencer",	"Type",
  "Sex",	"Age",	"Tissue",	"Extracted molecule",
  "Strain",	"Comments", "salmon", "transcriptome", "gtf")(ExtendedSample.apply)(ExtendedSample.unapply)



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

implicit val sampleCoder: HeaderCodec[Sample] = HeaderCodec.caseCodec(
  "GSM",	"GSE",	"Species",	"Sequencer",	"Type",
  "Sex",	"Age",	"Tissue",	"Extracted molecule",
  "Strain",	"Comments")(Sample.apply)(Sample.unapply)

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

object JSONReader extends JSONReader
trait JSONReader {

  def unsafeReadJson(references: Path): Json = {
    import io.circe.jackson.decode
    import io.circe.jackson.parse
    val str: String = read(references)
    val json = parse(str).toOption.get
    json
  }

  def unsafeReadField[A](json: Json, field: String)(implicit d: Decoder[A]) = json.hcursor.get[A](field).toOption.get

}

trait FolderExtractor {

  def extractNames(list: List[String]): List[String] = list.map(s=> Paths.get(s).getFileName.toString)
  def toFolder(folder: Path, value: String): String = if(value.trim=="" || value =="N/A") "" else (folder / value).toString()
  def toFolder(folder: Path, list: List[String]): List[String] = list.map(s=>toFolder(folder, s))

}


trait SalmonSample {
  def quant: String
  def expressions: String
  def compatible_fragment_ratio: String
  def exptected_lib_format: String
}

object FullSample extends FolderExtractor with JSONReader {

  def getLibRatio(path: Path): (String, Double) = {
    val json = unsafeReadJson(path)
    println()
    val ratio = unsafeReadField[Double](json, "compatible_fragment_ratio")
    println(ratio)
    val lib = unsafeReadField[String](json, "expected_format")
    println(lib)
    (lib, ratio)
  }

  lazy val sampleFilesHeaders = List(
    "forward_read_raw", "reverse_read_raw", "sra", //basic
    "forward_read_cleaned", "reverse_read_cleaned", "quality_html", "quality_json", //fastp
    "quant", "expressions", "compatible_fragment_ratio", "exptected_lib_format" //quantification
  )

  val headers =  ExtendedSample.headers ++sampleFilesHeaders

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

}

case class FullSample(gsm: String,	gse: String,	species: String,
                  sequencer: String, sample_type: String,	sex: String,
                  age: String,	tissue: String,	extracted_molecule: String,
                  strain: String,	comments: String, //Sample
                  salmon: String, transcriptome: String, gtf: String, //ExtendedSample
                  forward_read_raw: String, reverse_read_raw: String, sra: String,
                  forward_read_cleaned: String, reverse_read_cleaned: String, quality_html: String, quality_json: String,
                      quant: String, expressions: String, compatible_fragment_ratio: String, exptected_lib_format: String //, libType: String = "A"
                 ) extends PipelineSample with SalmonSample