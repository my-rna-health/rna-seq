task hello {
  String pattern
  File in

  command {
    echo 'hello ${pattern} world ${in}!'
  }

  output {
    #Array[String] matches = read_lines(stdout())
    String hello = "hello world"
  }
}

workflow wf {
  call hello
}