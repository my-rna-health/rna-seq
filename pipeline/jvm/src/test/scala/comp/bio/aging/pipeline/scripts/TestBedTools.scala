package comp.bio.aging.pipeline.scripts


object TestBedTools  {

  /**
    * How to run this script:
    * 1). Install Ammonite Shell http://www.lihaoyi.com/Ammonite/#Ammonite-Shell
    * On linux it is just one command:
    *  sudo curl -L -o /usr/local/bin/amm https://git.io/vMF2M && sudo chmod +x /usr/local/bin/amm && amm
    * 2) Run with the shell:
    *   amm benchmark.scala
    */
  import scala.concurrent.duration.FiniteDuration
  import ammonite.ops._
  import ammonite.ops.ImplicitWd._
  import ammonite.main.Router.main
  import ammonite.runtime.tools.time

  case class Benchmark[R](name: String, duration: FiniteDuration, result: R) {

    lazy val totalSecs: Long = duration.toSeconds
    lazy val hours: Long = totalSecs / 3600
    lazy val minutes: Long = (totalSecs % 3600) / 60
    lazy val seconds: Long = (totalSecs % 3600) % 60

    private def t(num: Long): String = if(num<10) s"0${num}" else num.toString

    lazy val message = s"${name} took: ${t(hours)} : ${t(minutes)} : ${t(seconds)}"
  }

  def writeBenchmark(benchmark: Benchmark[CommandResult]) = {
    import java.util._
    import java.text.SimpleDateFormat
    val timeStamp = new SimpleDateFormat("yyyy-MM-dd_HH:mm:ss").format(Calendar.getInstance().getTime())
    val folder =  pwd / (s"benchmark_${timeStamp}")
    mkdir! folder
    write.write(folder / "stdout.txt", benchmark.result.out.lines.mkString("\n"))
    write.write(folder / "stderr.txt", benchmark.result.err.lines.mkString("\n"))
    write.write(folder / "duration.txt", benchmark.message)
    println(s"benchmark results were successfully written to ${folder.toIO.getAbsolutePath}")
  }

  @main
  def intersect(file1Name: String, file2Name: String) = {
    val name = s"bedtools intersect -a, ${file1Name} -b ${file2Name}"
    println("starting the benchmark:\n"+name)
    val (result: CommandResult, duration: FiniteDuration) = time{
      %%("bedtools", "intersect", "-a", file1Name, "-b", file2Name)
    }
    writeBenchmark(Benchmark(name, duration, result))
  }

  @main
  def coverage(file1Name: String, file2Name:String, sorted: Boolean = false, counts: Boolean = false) = {
    val params = (if(sorted) " -sorted" else "") + (if(counts) " -counts" else "")
    val name = s"bedtools coverage -a, ${file1Name} -b ${file2Name}" + params
    println("starting the benchmark:\n"+name)
    val (result, duration) = time{(sorted, counts) match {
      case (true, true)  =>  %%("bedtools", "coverage", "-sorted", "-a", file1Name, "-b", file2Name, "-counts")
      case (true, false)  => %%("bedtools", "coverage", "-sorted", "-a", file1Name, "-b", file2Name)
      case (false, true)  => %%("bedtools", "coverage", "-a", file1Name, "-b", file2Name, "-counts")
      case (false, false) => %%("bedtools", "coverage", "-a", file1Name, "-b", file2Name)
    }}
    writeBenchmark(Benchmark(name, duration, result))
  }

}
