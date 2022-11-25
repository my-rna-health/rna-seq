version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

struct ExtractedRun {
    String run
    String folder
    Boolean is_paired
    Array[File] cleaned_reads
    Array[File] report
}

workflow extract_run{
    input {
        String layout
        String run
        String folder
        Boolean copy_cleaned = false
        Int extract_threads = 4
        Boolean aspera_download = true
        Boolean skip_technical = true
        Boolean original_name = false
    }
    Boolean is_paired = (layout != "SINGLE")

    call download { input: sra = run, aspera_download = aspera_download }
    call extract {input: sra = download.out, is_paired = is_paired, threads = extract_threads, skip_technical = skip_technical, original_name = original_name}
    call fastp { input: reads = extract.out, is_paired = is_paired }
    call files.copy as copy_report {
        input:
            destination = folder + "/report",
            files = [fastp.report_json, fastp.report_html]
    }
    if(copy_cleaned)
    {
        call files.copy as copy_cleaned_reads {
            input:
                destination = folder + "/reads",
                files = fastp.reads_cleaned
        }
    }

    output {
        ExtractedRun out = object {run: run, folder: folder, is_paired: is_paired, cleaned_reads: fastp.reads_cleaned, report: copy_report.out}
    }
}


task download {
    input {
        String sra
        Boolean aspera_download
    }
    #prefetch --ascp-path "/root/.aspera/connect/bin/ascp|/root/.aspera/connect/etc/asperaweb_id_dsa.openssh" --force yes -O results ~{sra}
    command {
        ~{if(aspera_download) then "download_sra_aspera.sh " else "prefetch -X 9999999999999 --force yes -O results -t http "} ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/comp-bio-aging/download_sra:latest"
        maxRetries: 1
    }

    output {
        File? a = "results" + "/" + sra + ".sra"
        File? b = "results" + "/" + sra + "/" + sra + ".sra"
        File out = select_first([a, b])
    }
}

task extract {
    input {
        File sra
        Boolean is_paired
        Boolean skip_technical
        Int threads
        Boolean original_name = false

    }

    String name = basename(sra, ".sra")
    String folder = "extracted"
    String prefix = folder + "/" + name

    #see https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump for docs

    command {
        ~{if(original_name) then "fastq_dump --origfmt " else "fasterq-dump "+ "--threads " + threads +" --progress "} --outdir ~{folder} --split-files ~{if(skip_technical) then "--skip-technical" else ""} ~{sra}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/download_sra:latest"
        maxRetries: 2
    }

    output {
        Array[File] out = glob(prefix+"*")
        #Array[File] out = if(is_paired) then [prefix + "_1.fastq",  prefix + "_2.fastq"] else [prefix + ".fastq"]
    }
}


task fastp {
    input {
        Array[File] reads
        Boolean is_paired
    }

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
        -i ~{reads[0]} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
        ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:2489fe56260bde05bdf72a8ead4892033b9a05dc4525affb909405bea7839d1b" #0.23.2--h5f740d0_3
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
                                    then [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz", basename(reads[1], ".fastq.gz") + "_cleaned.fastq.gz"]
                                    else [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz"]
    }
}