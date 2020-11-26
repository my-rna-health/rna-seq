version development

workflow rna_spades {

    input {
        Array[File] reads
        String destination
        Int max_memory = 55
        Int threads = 22
    }


    call fastp { input: reads = reads }

    call copy as copy_report {
         input:
            destination = destination + "/report",
            files = [fastp.report_json, fastp.report_html]
        }

    call rna_spades {
        input: reads = fastp.reads_cleaned, max_memory  = max_memory, threads = threads
    }

    call copy as copy_results {
         input:
            destination = destination,
            files = [rna_spades.out]
        }

    output {
        File out = rna_spades.out
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


task rna_spades {
    input {
        String results = "results"
        Array[File] reads
        String max_memory = 55
        Int threads
    }
    command {
         rnaspades.py -1 ~{reads[0]} -2 ~{reads[1]} --memory ~{max_memory} --threads ~{threads} -o ~{results}
    }

    runtime {
        docker: "quay.io/biocontainers/spades@sha256:9fc72d13bdd3b33af6c8f9bf03512dc486a50957d41eb27ed98eca0b98fa50ba"#:3.14.0--h2d02072_0"
        docker_memory: "${max_memory+1}G"
        docker_cpu: "${threads}"
    }

    output {
        File out = "results"
    }
}


task rna_quast {

    input {
        File contigs
        File? reference
        File? features
        Int? threads = 4
        File? features
        String output_folder = "results"
        Int min_contig = 50
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