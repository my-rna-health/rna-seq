workflow samtools {
    String name
    Array[File] aligments
    Int threads

     scatter (file in aligments) {
        call convert {
            input:
                sam = file,
                name = basename(file, ".sam"),
                threads = threads
        }
     }

     call merge {
        input:
            bams = convert.out,
            name = name
     }

     output {
        File merged = merge.out
     }

}

task convert {
    File sam
    String name
    Int threads

    command {
        samtools view -h -bS ${sam} | samtools sort - -@ ${threads} -o ${name}.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:b3bb39957750bc3c448e22488e75a7ec17fad03c1c9ab4a76aa2a44fc3843b36"
    }

    output {
        File out = name + ".bam"
    }
}

task merge {

    String name
    Array[File] bams

    command {
        samtools merge ${name}.bam ${sep=' ' bams}
    }

    runtime {
        docker: "biocontainers/samtools@sha256:b3bb39957750bc3c448e22488e75a7ec17fad03c1c9ab4a76aa2a44fc3843b36"
    }

    output {
        File out = name + ".bam"
    }
}