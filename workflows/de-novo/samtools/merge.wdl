workflow merge {
    String name
    Array[File] aligments

     scatter (file in aligments) {
        call convert {
            input: sam = file, name = basename(file)
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
    File name

    command {
        samtools view -bS file.sam | samtools sort - ${name}
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