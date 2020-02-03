version development

workflow idba {
    input {
        Array[File]+ reads
        String destination
        File reference
        Array[File] adapters
        Int q = 35
        Boolean short
    }

    call fastp { input: reads = reads,
        adapters = adapters, q = q
    }

    call idba_hybrid{
        input:
            reads = fastp.reads_cleaned, reference = reference, short = short
    }

    call copy{
        input: files = [idba_hybrid.out], destination = destination
    }


}



task fastp {
    input {
        Array[File] reads
        Array[String] adapters
        Int q = 35
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis --correction \
            --adapter_sequence ~{adapters[0]} --adapter_sequence_r2 ~{adapters[1]} -q ~{q} \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fq")}_cleaned.fq \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fq") +"_cleaned.fq" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:ac9027b8a8667e80cc1661899fb7e233143b6d1727d783541d6e0efffbb9594e" #0.20.0--hdbcaa40_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fq") + "_cleaned.fq", basename(reads[1], ".fq") + "_cleaned.fq"]
            else [basename(reads[0], ".fq") + "_cleaned.fq"]
    }
}

task idba_hybrid {
    input {
        Array[File] reads
        File reference
        Boolean short = true
    }

    command {
        fq2fa --merge ~{reads[0]} ~{reads[1]} reads_merged.fa
        idba_hybrid ~{if(short) then "--read" else "--long_read"} reads_merged.fa -o results --reference ~{reference}
    }

    runtime {
        docker: "quay.io/biocontainers/idba@sha256:51291ffeeecc6afab8d56bf33dffd0c2cb5e24d8a545a5ea93bb795d6af12fa0" #1.1.3--1
    }

    output {
        File out = "results"
    }
}



task fastp_merge {
    input {
        Array[File] reads
        Array[String] adapters
        Int q = 35
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis --correction \
            --adapter_sequence ~{adapters[0]} --adapter_sequence_r2 ~{adapters[1]} -q ~{q} \
            --merge --merged_out reads_merged.fq \
            -i ~{reads[0]} \
             ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1] else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:ac9027b8a8667e80cc1661899fb7e233143b6d1727d783541d6e0efffbb9594e" #0.20.0--hdbcaa40_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        File merged_reads = "reads_merged.fq"
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