#! /usr/local/bin/amm
@ interp.load.ivy("com.lihaoyi" %% "fastparse" % "0.4.3")

import ammonite.ops._
import ammonite.ops.ImplicitWd._

import fastparse.all._

@main
def main(file: Path, where: Path = pwd) = {
  val to = if(where.isFile) where else where / "result.fasta"
  val strings = read.lines! file
  val frag = "Longest kmer: "
  val lines = {
    for{
    s <- strings
    if s.contains(frag)
    txt = s.substring(s.indexOf(frag) + frag.length)
  } yield s">${frag+txt}\n${txt}\n"
  }
  write(to, lines)
  //println(s.substring(s.indexOf(frag)))
}
