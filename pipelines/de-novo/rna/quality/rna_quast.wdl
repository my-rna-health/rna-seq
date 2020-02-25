version development

workflow quast{
    input {
       File genome
       File gtf
       File transcripts
       Array[File] reads
       Int? threads
       String destination
       String output_name = "results"
    }

    call rna_quast {
        input:
            genome = genome, gtf = gtf, transcripts = transcripts, reads = reads, output_folder = output_name
    }


    call copy as copy_results {
         input:
            destination = destination,
            files = [rna_quast.out]
        }

    output {
        File out = copy_results.out[0]
    }

}

task rna_quast {

    input {
        File genome
        File gtf
        File transcripts
        Array[File] reads
        String output_folder = "results"
    }

    command {
        rnaQUAST.py --transcripts ~{transcripts} --reference ~{genome} --gtf ~{gtf} -1 ~{reads[0]} -2 ~{reads[1]}
    }

    runtime {
        docker: "quay.io/biocontainers/rnaquast@sha256:116b3f0ae54a18e10708dfe9a4edeec2f244b61a81f92b73e0d5392649b4dc75" #:2.0.0--0
    }

    output {
        File out = output_folder
    }
}

task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{where}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}