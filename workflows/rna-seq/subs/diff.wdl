workflow diff {

 Array[Pair[String, Array[File]]] conditions

  scatter (condition in conditions) {
      call condition  {
          input:
              condition = condition.left,
              samples =  condition.right
      }
  }
}

task trace {

    String condition
    Array[File] samples

    command {
     echo condition=${condition}
    }

    output {
       Array[String] out = read_lines(stdout())
    }
}