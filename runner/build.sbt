organization := "comp.bio.aging"

name := "rna-seq-runner"

scalaVersion :=  "2.12.3"

version := "0.0.1"

resolvers += sbt.Resolver.bintrayRepo("comp-bio-aging", "main")

resolvers += sbt.Resolver.bintrayRepo("denigma", "denigma-releases")

resolvers += "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/"

resolvers += "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/"

licenses += ("MPL-2.0", url("http://opensource.org/licenses/MPL-2.0"))

isSnapshot := true

scalacOptions ++= Seq( "-target:jvm-1.8", "-feature", "-language:_" )

javacOptions ++= Seq("-source", "1.8", "-target", "1.8", "-Xlint", "-J-Xss5M", "-encoding", "UTF-8")

libraryDependencies  ++= Seq(
  "comp.bio.aging" %% "cromwell-client" % "0.0.6.1"
)
