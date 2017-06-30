workflow test {

  File conditions

  Array[Array[String]] series =  read_tsv(conditions)

  scatter (row in series) {
    call trace {
        input: row = row
    }
  }
}

task trace {

    Array[String] row

    command {
     echo ${sep=" " row}
    }

    output {
       Array[String] out = read_lines(stdout())
    }
}



