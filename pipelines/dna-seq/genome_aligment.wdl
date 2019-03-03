version 1.0

workflow genome_alignment{
    input {

        Array[File] reads
        Array[String] adapters

        Int threads

        File reference
        String results_folder
    }

    call report as initial_report_1_call {
      input:
        sampleName = basename(reads[0], ".fastq.gz"),
        file = reads[0]
      }

    call report as initial_report_2_call {
      input:
        sampleName = basename(reads[1], ".fastq.gz"),
        file = reads[1]
      }

    call copy as copy_initial_quality_reports {
    input:
        files = [initial_report_1_call.out, initial_report_2_call.out],
        destination = results_folder + "/quality/initial/"
    }

    call atropos_illumina_trim as atropos_illumina_trim_call {
      input:
        reads = reads,
        adapters = adapters,
        threads = threads
    }

    call copy as copy_trimmed {
    input:
        files = [atropos_illumina_trim_call.out[0], atropos_illumina_trim_call.out[1]],
        destination = results_folder + "/trimmed/"
    }

    call report as final_report_1_call {
        input:
          sampleName = basename(atropos_illumina_trim_call.out[0], ".fastq.gz"),
          file = atropos_illumina_trim_call.out[0]
        }

    call report as final_report_2_call {
        input:
          sampleName = basename(atropos_illumina_trim_call.out[1], ".fastq.gz"),
          file = atropos_illumina_trim_call.out[1]
        }

    call copy as copy_cleaned_quality_reports {
    input:
        files = [final_report_1_call.out, final_report_2_call.out],
        destination = results_folder + "/quality/cleaned/"
    }


    call minimap2 {
        input:
            reads = atropos_illumina_trim_call.out,
            reference = reference
    }

    call samtools_conversion {
        input:
            sam = minimap2.out
    }

    call picard_readgroups_sort {
        input:
            bam = samtools_conversion.out
    }

    call picard_validation as picard_validation {
        input:
            bam = picard_readgroups_sort.out
    }

    call picard_indexbam as picard_indexbam {
        input:
            bam = picard_readgroups_sort.out
    }
}


task atropos_illumina_trim {
   input {
    Array[File] reads
    Array[String] adapters
    Int threads
    Int q = 18
    Float e = 0.1

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
    docker: "jdidion/atropos@sha256:c2018db3e8d42bf2ffdffc988eb8804c15527d509b11ea79ad9323e9743caac7"
  }

  output {
    Array[File] out = [basename(reads[0], ".fastq.gz") + "_trimmed.fastq.gz",  basename(reads[1], ".fastq.gz") + "_trimmed.fastq.gz"]
  }
}

task minimap2 {
    input {
        Array[File] reads
        File reference
    }

    command {
        minimap2 \
        -ax \
        sr \
        -L \
        ~{reference} \
        ~{reads[0]} \
        ~{reads[1]} \
        > aln.sam
    }

    runtime {
        docker: "genomicpariscentre/minimap2@sha256:536d7cc40209d4fd1b700ebec3ef9137ce1d9bc0948998c28b209a39a75458fa"
      }
    output {
      File out = "aln.sam"
    }
}

task samtools_conversion {
    input {
        File sam
    }

    command {
        samtools \
        view \
        -bS \
        ~{sam} \
        > aln.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:6644f6b3bb8893c1b10939406bb9f9cda58da368100d8c767037558142631cf3"
    }

    output {
        File out = "aln.bam"
    }

}

task picard_readgroups_sort{
    input {
        File bam
    }
    command {
picard AddOrReplaceReadGroups \
I=~{bam} \
O=aln2.bam \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=20 \
SORT_ORDER=coordinate
    }

    runtime {
        docker: "biocontainers/picard:v2.3.0_cv3"
    }

    output {
        File out = "aln2.bam"
    }

}

task picard_validation{
    input {
        File bam
    }
    command {
    picard ValidataeSamFile \
    I=~{bam} \
    O=log.txt \
    MODE=SUMMARY
    }

    runtime {
        docker: "biocontainers/picard:v2.3.0_cv3"
      }

    output {
        File out = "log.txt"
    }

}

task picard_indexbam {
    input {
        File bam
    }

    command {
        picard BuildBamIndex \
        INPUT=~{bam}
    }

    runtime {
        docker: "biocontainers/picard:v2.3.0_cv3"
      }

    output {
        File out = "aln2.bai"
      }

}

task report {
  input {
    String sampleName
    File file
  }

  command {
    /opt/FastQC/fastqc ~{file} -o .
  }

  runtime {
    docker: "quay.io/ucsc_cgl/fastqc@sha256:86d82e95a8e1bff48d95daf94ad1190d9c38283c8c5ad848b4a498f19ca94bfa"
  }

  output {
    File out = sampleName+"_fastqc.zip"
  }
}

task copy {
    input {
        Array[File] files
        String destination
    }

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{sep=' ' files} ~{destination}
    }

    output {
        Array[File] out = files
    }
}
