import coursier.maven.MavenRepository

import scala.util.{Failure, Try}
interp.repositories() ++= Seq(
  MavenRepository("https://dl.bintray.com/comp-bio-aging/main/")
)
@
import $ivy.`com.github.pathikrit::better-files:3.8.0`
import $ivy.`com.nrinaudo::kantan.csv:0.5.1`
import $ivy.`com.nrinaudo::kantan.csv-generic:0.5.1`
import $ivy.`group.aging-research::geo-fetch:0.0.6`
println("Script dependencies loaded!")