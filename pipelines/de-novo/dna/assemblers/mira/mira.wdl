version development

workflow mira {

    input {
        Array[File] reads
        String destination
        String project_name
        File? reference
        Int threads = 15
    }

    call fastp { input: reads = reads }

    call copy as copy_report {
         input:
            destination = destination + "/report",
            files = [fastp.report_json, fastp.report_html]
        }

    call mira_de_novo {
        input: reads = fastp.reads_cleaned, project_name = project_name, threads = threads
    }

    call copy as copy_results {
         input:
            destination = destination,
            files = [mira_de_novo.out]
        }

    output {
        File out = mira_de_novo.out
    }

}


task mira_de_novo {
    input {
        Array[File] reads
        String project_name
        Int threads
    }

    command {
       echo "project = ~{project_name}
       job = genome,denovo,accurate
       readgroup = DataIlluminaPairedLib
       autopairing
       data = fastq::~{reads[0]} fastq::~{reads[1]}
       technology = solexa
       " > manifest.txt
       mira -t ~{threads} manifest.txt
    }

    runtime {
        docker: "quay.io/comp-bio-aging/mira"
    }

    output {
        File out = 	project_name + "_assembly"
    }
}


task fastp {
    input {
        Array[File] reads
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fq")}_cleaned.fq \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fq") +"_cleaned.fq" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5" #0.20.1--h8b12597_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fq") + "_cleaned.fq", basename(reads[1], ".fq") + "_cleaned.fq"]
            else [basename(reads[0], ".fq") + "_cleaned.fq"]
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