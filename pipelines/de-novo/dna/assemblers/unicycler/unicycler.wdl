version development

workflow unicycler {

    input {
        Array[File] reads
        String destination
        String mode = "bold"
    }

    call fastp { input: reads = reads }

    call copy as copy_report {
         input:
            destination = destination + "/report",
            files = [fastp.report_json, fastp.report_html]
        }

    call unicycler {
        input: reads = reads, mode = mode
    }

    call copy as copy_results {
         input:
            destination = destination,
            files = [unicycler.out]
        }

    output {
        File out = unicycler.out
    }

}


task unicycler {
    input {
        Array[File] reads
        String mode #= "bold"
    }

    command {
        unicycler --mode ~{mode} --short1 ~{reads[0]} --short2 ~{reads[1]} --out results
    }

    runtime {
        docker: "quay.io/biocontainers/unicycler@sha256:134af9199ba08e177ad5e01a28e9ae160ceb2f35235da76d4d1dfb2c62c41f2f" #0.4.8--py37h8b12597_0
    }

    output {
        File out = "results"
    }
}


task fastp {
    input {
        Array[File] reads
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fq.gz")}_cleaned.fq.gz \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fq.gz") +"_cleaned.fq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5" #0.20.1--h8b12597_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz", basename(reads[1], ".fq.gz") + "_cleaned.fq.gz"]
            else [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz"]
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