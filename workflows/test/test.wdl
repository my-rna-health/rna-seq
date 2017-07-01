import "tracer.wdl" as tracer

workflow test {

  File conditions

  Array[Array[String]] series =  read_tsv(conditions)

  scatter (row in series) {

   call tracer.tracer { input:
      row = row
    }
  }

  output {
    Array[Array[Array[String]]] out = tracer.out
  }
}

