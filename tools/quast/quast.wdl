version development

workflow quast{
    input {
       Array[File]+ contigs
       File? reference
       Int? threads
       String destination
    }

    call quast {
        input:
            contigs = contigs,
            reference = reference,
            threads = threads
    }


    call copy as copy_results {
         input:
            destination = destination,
            files = [quast.out]
        }

    output {
        File out = copy_results.out[0]
    }

}

task quast {

    input {
        Array[File]+ contigs
        File? reference
        Int? threads = 4
        File? features
        String output_folder = "results"
        Float min_contig = 50
        String? type
    }

    command {
        quast.py ~{if defined(reference) then "--reference " + reference else ""} \
         ~{if defined(threads) then "--threads " + threads else ""} ~{contigs} \
         --output ~{output_folder} \
         ~{if defined(features) then "--features " + features + (if(defined(type)) then "--type " + type else "") else "" } \
         --min-contig ~{min_contig} \
         ~{sep=" " contigs}
    }

    runtime {
        docker: "quay.io/biocontainers/quast@sha256:8924f9a568deaa58118f36e47e333534ccb760dd51ed61f3fbd68fde9864c7c4"
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