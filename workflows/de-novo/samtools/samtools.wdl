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
        samtools view -bS ${sam} | samtools sort - -@ ${threads} -o ${name}.bam
    }

    runtime {
        docker: "quay.io/comp-bio-aging/samtools@sha256:2654f4f9d46d9a75b8c1d9e44cf91128f512b78a29e9f9fb15106bb278f03437"
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
        docker: "quay.io/comp-bio-aging/samtools@sha256:2654f4f9d46d9a75b8c1d9e44cf91128f512b78a29e9f9fb15106bb278f03437"
    }

    output {
        File out = name + ".bam"
    }
}