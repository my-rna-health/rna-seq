version development

workflow plasmid_spades {

    input {
        Array[File] reads
        String destination
    }


    call fastp { input: reads = reads }

    call copy as copy_report {
         input:
            destination = destination + "/report",
            files = [fastp.report_json, fastp.report_html]
        }

    call plasmid_spades {
        input: reads = reads,
    }

    call copy as copy_results {
         input:
            destination = destination,
            files = [plasmid_spades.out]
        }

    output {
        File out = plasmid_spades.out
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

task plasmid_spades {
    input {
        String results = "results"
        Array[File] reads
        String cut_off = "auto"
    }
    command {
        plasmidspades.py -1 ~{reads[0]} -2 ~{reads[1]} --cov-cutoff ~{cut_off} -o ~{results}
    }

    runtime {
        docker: "quay.io/biocontainers/spades@sha256:9fc72d13bdd3b33af6c8f9bf03512dc486a50957d41eb27ed98eca0b98fa50ba"#:3.14.0--h2d02072_0"
    }

    output {
        File out = "results"
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