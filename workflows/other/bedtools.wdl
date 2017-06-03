workflow bedtest {

    File file1
    File file2

    call intersect {
        input: a = file1,
            b = file2
    }
}

task intersect {

  File a
  File b

  #bbedtools intersect -a cpg.bed -b exons.bed

  command {
    bedtools intersect -a ${a} -b ${b} > intersection.txt
  }

  runtime {
    docker: "biocontainers/bedtools"
  }

  output {
    File intersection = "intersection.txt"
  }
}