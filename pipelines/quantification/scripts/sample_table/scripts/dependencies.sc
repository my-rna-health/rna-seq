import coursierapi.MavenRepository

interp.repositories.update(
  interp.repositories() ::: List(MavenRepository.of("https://dl.bintray.com/comp-bio-aging/main/"))
)
@
import $ivy.`com.github.pathikrit::better-files:3.8.0`
import $ivy.`com.nrinaudo::kantan.csv:0.6.0`
import $ivy.`com.nrinaudo::kantan.csv-generic:0.6.0`
import $ivy.`group.aging-research::geo-fetch:0.0.11`
println("Script dependencies loaded!")