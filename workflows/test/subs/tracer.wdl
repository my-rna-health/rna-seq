workflow tracer {
    Array[String] row

     scatter (sample in row) {
            if(sample != row[0]) {
                call trace {
                    input:
                        condition = row[0],
                        sample = sample
                }
            }
        }

    output {
        Array[Array[String]] out = select_all(trace.out)
    }
}

task trace {

    String condition
    String sample

    command {
     echo condition=${condition} experiment=${sample}
    }

    output {
       Array[String] out = read_lines(stdout())
    }
}