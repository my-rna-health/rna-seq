version 1.0

workflow methylation {
  input {
        Array[File] reads
        Array[String] adapters
        Int threads

        String results_folder
    }
}


task atropos_illumina_trim {
   input {
      Array[File] reads
    Array[String] adapters
    Int threads
    Int q = 18
    Int e = 0.1

   }

  command {
    atropos trim \
    -a ~{adapters[0]} \
    -A ~{adapters[1]} \
    -pe1 ~{reads[0]} \
    -pe2 ~{reads[1]} \
    -o ~{basename(reads[0], ".fastq.gz")}_trimmed.fastq.gz \
    -p ~{basename(reads[1], ".fastq.gz")}_trimmed.fastq.gz \
    --minimum-length 35 \
    --aligner insert \
    -q ~{q} \
    -e ~{e} \
    --threads ~{threads} \
    --correct-mismatches liberal
  }

  runtime {
    docker: "jdidion/atropos:latest"
  }

  output {
    Array[File] out = [basename(reads[0], ".fastq.gz") + "_trimmed.fastq.gz",  basename(reads[1], ".fastq.gz") + "_trimmed.fastq.gz"]
  }
}