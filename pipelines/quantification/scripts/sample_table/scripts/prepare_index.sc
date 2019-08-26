#! /usr/local/bin/amm
import $exec.dependencies

import ammonite.ops._
import ammonite.ops.ImplicitWd._
import better.files._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._
import io.circe.optics.JsonPath.root
import io.circe._, io.circe.parser._
import scala.util._

import io.circe.optics.JsonPath._
/*
val title = root.title.string
val organism = root.organism.name.string
val taxid = root.organism.taxid.string
val sequencer = root.sequencer.string
val strategy = root.library.strategy.string
val selection = root.library.selection.string
val source = root.library.source.string
val extraction_source   = root.extraction.source.string
val molecule = root.extraction.molecule.string
val protocol = root.extraction.protocol.string
*/


