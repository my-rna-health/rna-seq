task download {

  File in

  command {
    /opt/sratoolkit
  }

  runtime {
    docker: "itsjeffreyy/sratoolkit"
  }

  output {
    #Array[String] matches = read_lines(stdout())
    "/opt/sratoolkit"
  }
}

workflow seq {
  call download
}