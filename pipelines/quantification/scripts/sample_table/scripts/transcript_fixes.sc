#! /usr/local/bin/amm
import $exec.dependencies
import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._


@main
def main(fasta: Path = Path("/data/ensembl/97/species/homo_sapiens/test_rna.txt"),
         output: Path = Path("/data/ensembl/97/species/homo_sapiens/test_rna_stable.txt")) = {
  println(s"deleting dots from ${fasta} and writing results to ${output}")
  val o = output.toIO.toScala.createIfNotExists().clear()
  fasta.toIO.toScala.lineIterator
  for{
    l <- fasta.toIO.toScala.lineIterator
  } {
    val n  =if(l.startsWith(">") && l.contains(".")) {
      l.indexOf(".")
      match {
        case -1 => l
        case i if i <= l.length -2 => l.substring(0, i) + l.substring(Math.max(i+2, l.indexOf(" ", i)))
        case _ => l
      }

    } else l
    o.appendLine(n)
  }
  println(s"fixing dot in ${fasta} finished, results saved at ${o.pathAsString}")
}