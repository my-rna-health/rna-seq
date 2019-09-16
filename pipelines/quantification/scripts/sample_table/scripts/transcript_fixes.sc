#! /usr/local/bin/amm
import $exec.dependencies
import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._


@main
def main(fasta: Path = Path("/data/ensembl/97/species/homo_sapiens/test_rna.txt"), output: Path = Path("/")) = {
  //val lines = read.lines! fasta
  //lines.map(str => str.substring(0, ))
  //write(output, test_rna.txt
  val o = output.toIO.toScala.createIfNotExists()
  fasta.toIO.toScala.lineIterator
  for{
    l <- fasta.toIO.toScala.lineIterator
  } {
    val n  =if(l.startsWith(">") && l.contains(".")) {
      val i = l.indexOf(".")
      l.substring(0, i)
    } else l
    println(n)

  }
}