version development

struct NamedRun{
    String name
    String run
}

workflow download_and_merge{
    input {
        Array[NamedRun] named_runs
        String destination
        Boolean aspera_download = true
        Int extract_threads = 4
        Boolean is_paired = true
    }

    scatter (named_run in named_runs) {
        String name = named_run.name
        String run = named_run.run
        call download { input: sra = run, aspera_download = aspera_download }
        call extract {input: sra = download.out, is_paired = is_paired, threads = extract_threads}
        String folder = destination + "/" + name + "/" + run
        call copy as copy_forward {
                 input:
                    destination = folder + "/reads",
                    files = [extract.out[0]]
            }
        call copy as copy_reverse {
                 input:
                    destination = folder + "/reads",
                    files = [extract.out[1]]
            }

        call fastp_uncompressed { input: reads = extract.out, is_paired = is_paired }
        call copy as copy_report {
         input:
            destination = folder + "/report",
            files = [fastp_uncompressed.report_json, fastp_uncompressed.report_html]
        }
        call copy as copy_cleaned_reads {
             input:
                destination = folder + "/cleaned_reads",
                files = fastp_uncompressed.reads_cleaned
        }
      }

      call merge_reads {
        input:
            forward = flatten(copy_forward.out),
            reverse = flatten(copy_reverse.out)
      }

      call copy as copy_merged {
        input: files = merge_reads.out,
        destination = destination + "/" + "merged"
      }

      output {
        Array[File] merged = copy_merged.out
        Array[Array[File]] cleaned_reads = copy_cleaned_reads.out
      }

}



task download {
    input {
        String sra
        Boolean aspera_download
    }
    #prefetch --ascp-path "/root/.aspera/connect/bin/ascp|/root/.aspera/connect/etc/asperaweb_id_dsa.openssh" --force yes -O results ~{sra}
    command {
        ~{if(aspera_download) then "download_sra_aspera.sh " else "prefetch --force yes -O results -t http "} ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/comp-bio-aging/download_sra:master"
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
        Int threads
    }

    String name = basename(sra, ".sra")
    String folder = "extracted"
    String prefix = folder + "/" + name
    String prefix_sra = prefix + ".sra"

    #see https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump for docs

    command {
        fasterq-dump --outdir ~{folder} --threads ~{threads} --progress --split-files --skip-technical ~{sra}
        ~{if(is_paired) then "mv" + " " + prefix_sra + "_1.fastq" + " " + prefix + "_1.fastq"  else "mv" + " " + prefix_sra + ".fastq" + " " + prefix + ".fastq"}
        ~{if(is_paired) then "mv" + " " + prefix_sra + "_2.fastq" + " " + prefix + "_2.fastq"  else ""}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/download_sra:master"
        maxRetries: 1
    }

    output {
        Array[File] out = if(is_paired) then [prefix + "_1.fastq",  prefix + "_2.fastq"] else [prefix + ".fastq"]
     }
}

task fastp_uncompressed {
    input {
        Array[File] reads
        Boolean is_paired
    }

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fastq")}_cleaned.fastq \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq") +"_cleaned.fastq" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:ac9027b8a8667e80cc1661899fb7e233143b6d1727d783541d6e0efffbb9594e" #0.20.0--hdbcaa40_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fastq") + "_cleaned.fastq", basename(reads[1], ".fastq") + "_cleaned.fastq"]
            else [basename(reads[0], ".fastq") + "_cleaned.fastq"]
    }
}


task merge_reads {
    input {
        Array[File] forward
        Array[File] reverse
    }

    command {
        cat ~{sep=" " forward} > forward.fastq
        cat ~{sep=" " reverse} > reverse.fastq
    }

    output {
        Array[File] out = ["forward.fastq", "reverse.fastq"]
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